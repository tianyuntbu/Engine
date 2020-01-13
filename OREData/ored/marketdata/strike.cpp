/*
 Copyright (C) 2019 Quaternion Risk Management Ltd
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

#include <ored/marketdata/strike.hpp>
#include <ored/utilities/parsers.hpp>
#include <ored/utilities/to_string.hpp>
#include <boost/algorithm/string.hpp>

using namespace QuantLib;
using std::ostream;
using std::ostringstream;
using std::string;
using std::vector;

namespace ore {
namespace data {

AbsoluteStrike::AbsoluteStrike() {}

AbsoluteStrike::AbsoluteStrike(Real strike) : strike_(strike) {}

Real AbsoluteStrike::strike() const {
    return strike_;
}

void AbsoluteStrike::fromString(const string& strStrike) {
    strike_ = parseReal(strStrike);
}

string AbsoluteStrike::toString() const {
    return to_string(strike_);
}

DeltaStrike::DeltaStrike() : deltaType_(DeltaVolQuote::Spot), optionType_(Option::Call) {}

DeltaStrike::DeltaStrike(DeltaVolQuote::DeltaType deltaType, Option::Type optionType, Real delta)
    : deltaType_(deltaType), optionType_(optionType), delta_(delta) {}

DeltaVolQuote::DeltaType DeltaStrike::deltaType() const {
    return deltaType_;
}

Option::Type DeltaStrike::optionType() const {
    return optionType_;
}

Real DeltaStrike::delta() const {
    return delta_;
}

void DeltaStrike::fromString(const string& strStrike) {
    
    // Expect strStrike of form: DEL / Spot|Fwd|PaSpot|PaFwd / Call|Put / DELTA_VALUE
    vector<string> tokens;
    boost::split(tokens, strStrike, boost::is_any_of("/"));
    QL_REQUIRE(tokens.size() == 4, "DeltaStrike::fromString expects 4 tokens.");
    QL_REQUIRE(tokens[0] == "DEL", "DeltaStrike::fromString expects 1st token to equal 'DEL'.");
    
    deltaType_ = parseDeltaType(tokens[1]);
    optionType_ = parseOptionType(tokens[2]);
    delta_ = parseReal(tokens[3]);
}

string DeltaStrike::toString() const {

    // Write string of form: DEL / Spot|Fwd|PaSpot|PaFwd / Call|Put / DELTA_VALUE
    ostringstream oss;
    oss << "DEL/" << deltaType_ << "/" << optionType_ << "/" << to_string(delta_);
    return oss.str();
}

AtmStrike::AtmStrike() : atmType_(DeltaVolQuote::AtmSpot) {}

AtmStrike::AtmStrike(DeltaVolQuote::AtmType atmType, boost::optional<DeltaVolQuote::DeltaType> deltaType)
    : atmType_(atmType), deltaType_(deltaType) {
    check();
}

DeltaVolQuote::AtmType AtmStrike::atmType() const {
    return atmType_;
}

boost::optional<DeltaVolQuote::DeltaType> AtmStrike::deltaType() const {
    return deltaType_;
}

void AtmStrike::fromString(const string& strStrike) {

    // Expect strStrike of form:
    // "ATM / AtmSpot|AtmFwd|AtmDeltaNeutral|AtmVegaMax|AtmGammaMax|AtmPutCall50"
    // with an optional following "/ DEL / Spot|Fwd|PaSpot|PaFwd"
    vector<string> tokens;
    boost::split(tokens, strStrike, boost::is_any_of("/"));
    QL_REQUIRE(tokens.size() == 2 || tokens.size() == 4, "AtmStrike::fromString expects 2 or 4 tokens.");
    QL_REQUIRE(tokens[0] == "ATM", "AtmStrike::fromString expects 1st token to equal 'ATM'.");

    atmType_ = parseAtmType(tokens[1]);

    deltaType_ = boost::none;
    if (tokens.size() == 4) {
        QL_REQUIRE(tokens[2] == "DEL", "AtmStrike::fromString expects 3rd token to equal 'DEL'.");
        deltaType_ = parseDeltaType(tokens[3]);
    }

    check();
}

string AtmStrike::toString() const {
    
    // Write strike of form:
    // "ATM / AtmSpot|AtmFwd|AtmDeltaNeutral|AtmVegaMax|AtmGammaMax|AtmPutCall50"
    // with an optional following "/ DEL / Spot|Fwd|PaSpot|PaFwd" if delta type is populated.
    ostringstream oss;
    oss << "ATM/" << atmType_;
    
    if (deltaType_) {
        oss << "/DEL/" << (*deltaType_);
    }

    return oss.str();
}

void AtmStrike::check() const {
    QL_REQUIRE(atmType_ != DeltaVolQuote::AtmNull, "AtmStrike type must not be AtmNull.");
    if (atmType_ == DeltaVolQuote::AtmDeltaNeutral) {
        QL_REQUIRE(deltaType_, "If AtmStrike type is AtmDeltaNeutral, we need a delta type.");
    } else {
        QL_REQUIRE(!deltaType_, "If AtmStrike type is not AtmDeltaNeutral, delta type should not be given.");
    }
    if (atmType_ == DeltaVolQuote::AtmPutCall50) {
        QL_REQUIRE(deltaType_ && *deltaType_ == DeltaVolQuote::Fwd, 
            "If AtmStrike type is AtmPutCall50, delta type must be AtmFwd.");
    }
}

MoneynessStrike::MoneynessStrike() : type_(Type::Spot) {}

MoneynessStrike::MoneynessStrike(Type type, Real moneyness) : type_(type), moneyness_(moneyness) {}

MoneynessStrike::Type MoneynessStrike::type() const {
    return type_;
}

Real MoneynessStrike::moneyness() const {
    return moneyness_;
}

void MoneynessStrike::fromString(const string& strStrike) {

    // Expect strStrike of form "MNY / Spot|Fwd / MONEYNESS_VALUE"
    vector<string> tokens;
    boost::split(tokens, strStrike, boost::is_any_of("/"));
    QL_REQUIRE(tokens.size() == 3, "MoneynessStrike::fromString expects 3 tokens.");
    QL_REQUIRE(tokens[0] == "MNY", "MoneynessStrike::fromString expects 1st token to equal 'MNY'.");

    type_ = parseMoneynessType(tokens[1]);
    moneyness_ = parseReal(tokens[2]);
}

string MoneynessStrike::toString() const {

    // Write strike of form "MNY / Spot|Fwd / MONEYNESS_VALUE"
    ostringstream oss;
    oss << "MNY/" << type_ << "/" << to_string(moneyness_);
    return oss.str();
}

ostream& operator<<(ostream& os, DeltaVolQuote::DeltaType type) {
    switch (type) {
    case DeltaVolQuote::Spot:
        return os << "Spot";
    case DeltaVolQuote::Fwd:
        return os << "Fwd";
    case DeltaVolQuote::PaSpot:
        return os << "PaSpot";
    case DeltaVolQuote::PaFwd:
        return os << "PaFwd";
    default:
        QL_FAIL("Unknown delta type");
    }
}

ostream& operator<<(ostream& os, DeltaVolQuote::AtmType type) {
    switch (type) {
    case DeltaVolQuote::AtmNull:
        return os << "AtmNull";
    case DeltaVolQuote::AtmSpot:
        return os << "AtmSpot";
    case DeltaVolQuote::AtmFwd:
        return os << "AtmFwd";
    case DeltaVolQuote::AtmDeltaNeutral:
        return os << "AtmDeltaNeutral";
    case DeltaVolQuote::AtmVegaMax:
        return os << "AtmVegaMax";
    case DeltaVolQuote::AtmGammaMax:
        return os << "AtmGammaMax";
    case DeltaVolQuote::AtmPutCall50:
        return os << "AtmPutCall50";
    default:
        QL_FAIL("Unknown atm type");
    }
}

ostream& operator<<(ostream& os, MoneynessStrike::Type type) {
    switch (type) {
    case MoneynessStrike::Type::Spot:
        return os << "Spot";
    case MoneynessStrike::Type::Forward:
        return os << "Fwd";
    default:
        QL_FAIL("Unknown moneyness type");
    }
}

MoneynessStrike::Type parseMoneynessType(const string& type) {
    if (type == "Spot") {
        return MoneynessStrike::Type::Spot;
    } else if (type == "Fwd") {
        return MoneynessStrike::Type::Forward;
    } else {
        QL_FAIL("Moneyness type '" << type << "' not recognized");
    }
}

boost::shared_ptr<BaseStrike> parseStrike(const string& strStrike) {
    
    boost::shared_ptr<BaseStrike> strike;

    // Expect strStrike to either:
    // 1. have a single token which means that we have an absolute strike, or
    // 2. have multiple tokens beginning with one of DEL, ATM or MNY
    vector<string> tokens;
    boost::split(tokens, strStrike, boost::is_any_of("/"));

    if (tokens.size() == 1) {
        strike = boost::make_shared<AbsoluteStrike>();
    } else if (tokens[0] == "DEL") {
        strike = boost::make_shared<DeltaStrike>();
    } else if (tokens[0] == "ATM") {
        strike = boost::make_shared<AtmStrike>();
    } else if (tokens[0] == "MNY") {
        strike = boost::make_shared<MoneynessStrike>();
    } else {
        QL_FAIL("Could not parse strike string '" << strStrike << "'.");
    }

    strike->fromString(strStrike);

    return strike;
}

}
}
