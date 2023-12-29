#include <qle/models/hwparametrization.hpp>
#include <qle/models/hwconstantparametrization.hpp>
#include <qle/pricingengines/analyticnfhwswaption.hpp>

#include <boost/test/unit_test.hpp>
//#include <oret/toplevelfixture.hpp>
#include <test/toplevelfixture.hpp>

#include <algorithm>
#include <boost/assign/std/vector.hpp>

#include <qle/models/gaussian1dcrossassetadaptor.hpp>
#include <qle/models/irlgm1fpiecewiseconstanthullwhiteadaptor.hpp>
#include <qle/models/lgm.hpp>
#include <qle/pricingengines/numericlgmmultilegoptionengine.hpp>

#include <qle/pricingengines/mclgmswaptionengine.hpp>

#include <ql/currencies/europe.hpp>
#include <ql/indexes/ibor/euribor.hpp>
#include <ql/models/shortrate/onefactormodels/gsr.hpp>
#include <ql/pricingengines/swap/discountingswapengine.hpp>
#include <ql/pricingengines/swaption/gaussian1dswaptionengine.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/time/all.hpp>

#include <boost/timer/timer.hpp>

using namespace QuantLib;
using namespace QuantExt;
using namespace boost::unit_test_framework;
using namespace boost::assign;
using namespace std;


BOOST_AUTO_TEST_CASE(testHwSwaptionEngines) {

    BOOST_TEST_MESSAGE("Testing pricing of the Swaption Engine measure in Hull White model...");

    Calendar cal = TARGET();
    Date evalDate(5, February, 2016);
    Date effectiveDate(cal.advance(evalDate,2 * Days));
    Date startDate(cal.advance(effectiveDate, 1 * Years));
    Date maturityDate(cal.advance(startDate, 9 * Years));
    
    Settings::instance().evaluationDate() = evalDate;

    // Set up the yield curve
    Real initial_rate = 0.02;
    Handle<YieldTermStructure> yts(boost::make_shared<FlatForward>(evalDate, initial_rate, Actual365Fixed()));

    // Create a fixed-rate and floating-rate leg for a vanilla interest rate swap (IRS)
    Real nominal = 1.0;
    Rate fixedRate = 0.02;
    DayCounter fixedDayCount = Thirty360(Thirty360::BondBasis);
    Schedule fixedSchedule(startDate, maturityDate, 1 * Years, cal, ModifiedFollowing, ModifiedFollowing,
                           DateGeneration::Forward, false);

    Handle<YieldTermStructure> forwardCurve(yts); // Same yield curve for the forward rates
    DayCounter floatingDayCount = Actual360();
    Schedule floatingSchedule(startDate, maturityDate, 6 * Months, cal, ModifiedFollowing, ModifiedFollowing,
                              DateGeneration::Forward, false);

    /* sigma, a m x n matrix */
    Matrix sigma(1, 1, 0.01);
    Array kappa(1, 0.03);
    boost::shared_ptr<IrHwParametrization> HwParam = boost::make_shared<IrHwConstantParametrization>(EURCurrency(), yts, sigma, kappa);
    boost::shared_ptr<HwModel> hwm = boost::make_shared<HwModel>(HwParam);
    
    boost::shared_ptr<PricingEngine> engine_nfhw = boost::make_share<AnalyticNfhwSwaptionEngine>(hwm, yts);
    //set the underlying swap -> swaption -> analytic engine
    boost::shared_ptr<IborIndex> euribor6m(boost::make_shared<Euribor>(6 * Months, yts));
    Spread spread = 0.0;
    boost::shared_ptr<VanillaSwap> undlSwap = boost::make_shared<VanillaSwap>(
            VanillaSwap(VanillaSwap::Payer, nominal, fixedSchedule, fixedRate, fixedDayCount, floatingSchedule, euribor6m,
                        spread, floatingDayCount));
    boost::shared_ptr<Swaption> swaption_nfhw = boost::make_shared<Swaption>(undlSwap, exercise);
    swaption_nfhw->setPricingEngine(engine_nfhw);
    Real nfhw_swaption_value = swaption_nfhw->NPV();
    
    boost::shared_ptr<StochasticProcess> process = hwm->stateProcess();

    // path generation
    Size n = 1000;   // number of paths
    Size seed = 121; // seed
    // maturity of test payoffs
    Time T = 1.0;
    // take large steps, but not only one (since we are testing)
    Size steps = static_cast<Size>(T * 24);
    TimeGrid grid(T, steps);
    PseudoRandom::rsg_type sg = PseudoRandom::make_sequence_generator(steps, seed);
    PathGenerator<PseudoRandom::rsg_type> pg(process, grid, sg, false);

    std::vector<Real> swaption_arr;
    
    for (Size i = 0; i < n; ++i) {
        Sample<Path> path = pg.next();
        Real fixed_leg = 0.0;
        for (Size j = 0; j < 9; ++j) {
            Date paymentDate(cal.advance(fixedSchedule[j]));
            Real fixed_df = discountBond(startDate, paymentDate, path.value, yts);
            Real fixed_pv = fixedRate * fixed_df;
            fixed_leg += fixed_pv;
        }
        Real float_leg = 0.0;
        Date d1 = startDate;
        for (Size k = 0; k < 18; ++k){
            Date d2(cal.advance(floatingSchedule[k]));
            DiscountFactor disc1 = yts->discount(d1);
            DiscountFactor disc2 = yts->discount(d2);
            Real floateRate = (disc1 / disc2 - 1.0) / 0.5;
            Real floate_df = discountBond(startDate, d2, path.value, yts);
            Real float_pv = floateRate * 0.5 * floate_df;
            float_leg += float_pv;
            d1 = d2;
        }
        Real swap_receiver_value = max(fixed_leg - float_leg,0);
        Real swaption_numeraire = numeraire(startDate, path.value, yts);
        Real swaption_value = swap_receiver_value / swaption_numeraire;
        swaption_arr.push_back(swaption_value);
    }
    Real sum = std::accumulate(swaption_arr.begin(), swaption_arr.end(), 0.0);
    Real mean = sum / swaption_arr.size();
    Real mc_swaption_value = mean;
    
    // tolerance for comparison analytic nfhw swaption engine vs unit test mc engine
    Real tol1 = 1.0E-3;
    if (std::abs(nfhw_swaption_value - mc_swaption_value) > tol) {
        BOOST_ERROR("inconsistent swaption npvs (analytic_nfhw="
                    << nfhw_swaption_value << ", mc_nfhw=" << mc_swaption_value;
    }
}

