/*
 Copyright (C) 2020 Quaternion Risk Management Ltd
 All rights reserved.
*/

// clang-format off
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
// clang-format on

#include <ored/portfolio/cbo.hpp>
#include <ored/portfolio/builders/cbo.hpp>
#include <oret/toplevelfixture.hpp>

#include <oret/datapaths.hpp>
#include <ored/configuration/curveconfigurations.hpp>
#include <ored/marketdata/csvloader.hpp>
#include <ored/marketdata/todaysmarket.hpp>

#include <ored/portfolio/portfolio.hpp>

#include <iostream>
#include <iomanip>

using namespace ore::data;
using namespace QuantExt;

BOOST_FIXTURE_TEST_SUITE(OREPlusCreditTestSuite, ore::test::TopLevelFixture)

BOOST_AUTO_TEST_SUITE(CBOTest)

BOOST_AUTO_TEST_CASE(testSimpleCBO) {
    BOOST_TEST_MESSAGE("Testing simple CBO...");

    ORE_REGISTER_REFERENCE_DATUM("CBO", CboReferenceDatum, true)
    ORE_REGISTER_ENGINE_BUILDER(CboMCEngineBuilder, false)
    ORE_REGISTER_TRADE_BUILDER("CBO", ore::data::CBO, false)

    Settings::instance().evaluationDate() = Date(31, Dec, 2018);
    Date asof = Settings::instance().evaluationDate();

    // Market
    auto conventions = boost::make_shared<Conventions>();
    conventions->fromFile(TEST_INPUT_FILE("conventions.xml"));
    InstrumentConventions::instance().setConventions(conventions);

    auto todaysMarketParams = boost::make_shared<TodaysMarketParameters>();
    todaysMarketParams->fromFile(TEST_INPUT_FILE("todaysmarket.xml"));
    auto curveConfigs = boost::make_shared<CurveConfigurations>();
    curveConfigs->fromFile(TEST_INPUT_FILE("curveconfig.xml"));
    auto loader = boost::make_shared<CSVLoader>(TEST_INPUT_FILE("market.txt"), TEST_INPUT_FILE("fixings.txt"), false);
    auto market = boost::make_shared<TodaysMarket>(asof, todaysMarketParams, loader, curveConfigs, false);

    // Portfolio to test market
    boost::shared_ptr<EngineData> engineData = boost::make_shared<EngineData>();
    engineData->fromFile(TEST_INPUT_FILE("pricingengine.xml"));
    boost::shared_ptr<EngineFactory> factory = boost::make_shared<EngineFactory>(engineData, market);

    Portfolio p;
    BOOST_CHECK_NO_THROW(p.fromFile(TEST_INPUT_FILE("cbo.xml")));
    BOOST_CHECK_NO_THROW(p.build(factory));

    // Pricing comparison
    double expectedNpv = 3013120.939;
    const Real tol = 0.01;

    BOOST_CHECK_NO_THROW(p.get("CBO-Constellation")->instrument()->NPV());
    BOOST_TEST_MESSAGE(p.get("CBO-Constellation")->instrument()->NPV());
    BOOST_CHECK_CLOSE(p.get("CBO-Constellation")->instrument()->NPV(), expectedNpv, tol);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
