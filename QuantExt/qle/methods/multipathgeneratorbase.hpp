/*
 Copyright (C) 2016 Quaternion Risk Management Ltd
 All rights reserved.

 This file is part of ORE, a free-software/open-source library
 for transparent pricing and risk analysis - http://opensourcerisk.org

 ORE is free software: you can redistribute it and/or modify it
 under the terms of the Modified BSD License.  You should have received a
 copy of the license along with this program.
 The license is also available online at <http://opensourcerisk.org>

 This program is distributed on the basis that it will form a useful
 contribution to risk analytics and model standardisation, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE. See the license for more details.
*/

/*! \file multipathgeneratorbase.hpp
    \brief base class for multi path generators
    \ingroup methods
*/

#pragma once

#include <ql/math/randomnumbers/rngtraits.hpp>
#include <ql/methods/montecarlo/brownianbridge.hpp>
#include <ql/methods/montecarlo/multipath.hpp>
#include <ql/methods/montecarlo/multipathgenerator.hpp>
#include <ql/methods/montecarlo/sample.hpp>
#include <ql/models/marketmodels/browniangenerators/sobolbrowniangenerator.hpp>
#include <ql/stochasticprocess.hpp>

namespace QuantExt {
using namespace QuantLib;

enum SequenceType { MersenneTwister, MersenneTwisterAntithetic, Sobol, SobolBrownianBridge };

//! Multi Path Generator Base
/*! \ingroup methods
 */
class MultiPathGeneratorBase {
public:
    virtual ~MultiPathGeneratorBase() {}
    virtual const Sample<MultiPath>& next() const = 0;
    virtual void reset() = 0;
};

//! Instantiation of MultiPathGenerator with standard PseudoRandom traits
/*! \ingroup methods
 */
class MultiPathGeneratorMersenneTwister : public MultiPathGeneratorBase {
public:
    MultiPathGeneratorMersenneTwister(const boost::shared_ptr<StochasticProcess>&, const TimeGrid&, BigNatural seed = 0,
                                      bool antitheticSampling = false);
    const Sample<MultiPath>& next() const override;
    void reset() override;

private:
    const boost::shared_ptr<StochasticProcess> process_;
    TimeGrid grid_;
    BigNatural seed_;

    boost::shared_ptr<MultiPathGenerator<PseudoRandom::rsg_type>> pg_;
    boost::shared_ptr<PathGenerator<PseudoRandom::rsg_type>> pg1D_;
    bool antitheticSampling_;
    mutable bool antitheticVariate_;
    mutable Sample<MultiPath> next_;
};

class MultiPathGeneratorMersenneTwisterAntithetic : public MultiPathGeneratorMersenneTwister {
public:
    MultiPathGeneratorMersenneTwisterAntithetic(const boost::shared_ptr<StochasticProcess>& p, const TimeGrid& grid,
                                                BigNatural seed = 0)
        : MultiPathGeneratorMersenneTwister(p, grid, seed, true) {}
};

//! Instantiation of MultiPathGenerator with standard LowDiscrepancy traits
/*! no Brownian bridge provided, use MultiPathGeneratorSobolBrownianBridge for this,
for the use of the seed, see ql/math/randomnumbers/sobolrsg.cpp

    \ingroup methods
*/
class MultiPathGeneratorSobol : public MultiPathGeneratorBase {
public:
    MultiPathGeneratorSobol(const boost::shared_ptr<StochasticProcess>&, const TimeGrid&, BigNatural seed = 0,
                            SobolRsg::DirectionIntegers directionIntegers = SobolRsg::JoeKuoD7);
    const Sample<MultiPath>& next() const override;
    void reset() override;

private:
    const boost::shared_ptr<StochasticProcess> process_;
    TimeGrid grid_;
    BigNatural seed_;
    SobolRsg::DirectionIntegers directionIntegers_;

    boost::shared_ptr<MultiPathGenerator<LowDiscrepancy::rsg_type>> pg_;
    boost::shared_ptr<PathGenerator<LowDiscrepancy::rsg_type>> pg1D_;
    mutable Sample<MultiPath> next_;
};

//! Instantiation using SobolBrownianGenerator from  models/marketmodels/browniangenerators
/*! \ingroup methods
 */
class MultiPathGeneratorSobolBrownianBridge : public MultiPathGeneratorBase {
public:
    MultiPathGeneratorSobolBrownianBridge(const boost::shared_ptr<StochasticProcess>&, const TimeGrid&,
                                          SobolBrownianGenerator::Ordering ordering = SobolBrownianGenerator::Steps,
                                          BigNatural seed = 0,
                                          SobolRsg::DirectionIntegers directionIntegers = SobolRsg::JoeKuoD7);
    const Sample<MultiPath>& next() const override;
    void reset() override;

private:
    const boost::shared_ptr<StochasticProcess> process_;
    TimeGrid grid_;
    SobolBrownianGenerator::Ordering ordering_;
    BigNatural seed_;
    SobolRsg::DirectionIntegers directionIntegers_;
    boost::shared_ptr<SobolBrownianGenerator> gen_;
    mutable Sample<MultiPath> next_;
    boost::shared_ptr<StochasticProcess1D> process1D_;
};

//! Make function for path generators
boost::shared_ptr<MultiPathGeneratorBase>
makeMultiPathGenerator(const SequenceType s, const boost::shared_ptr<StochasticProcess>& process,
                       const TimeGrid& timeGrid, const BigNatural seed,
                       const SobolBrownianGenerator::Ordering ordering = SobolBrownianGenerator::Steps,
                       const SobolRsg::DirectionIntegers directionIntegers = SobolRsg::JoeKuoD7);

//! Output function
std::ostream& operator<<(std::ostream& out, const SequenceType s);

} // namespace QuantExt
