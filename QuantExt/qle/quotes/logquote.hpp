/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2012 - 2016 Quaternion Risk Management Ltd.
 All rights reserved
*/

/*! \file logquote.hpp
    \brief stores log of quote for log-linear interpolation 
*/

#ifndef quantext_logquote_hpp
#define quantext_logquote_hpp

#include <ql/quote.hpp>

using namespace QuantLib;

namespace QuantExt {

    //! Class for storing logs of quotes for log-linear interpolation.
    /*! \test the correctness of the returned values is tested by
            checking them against the log of the returned values of q_
    */
    class LogQuote : public Quote, public Observer {
      public:
        LogQuote(Handle<Quote> q) : q_(q) {
            registerWith(q);
            update();
        }
        //! \name Inspectors
        //@{
        Real quote() const {
            return q_;
        }
        //@}
        //! \name Quote interface
        //@{
        Real value() const {
            return logValue_;
        }
        bool isValid() const {
            return q_->isValid();
        }
        //@}
        //! \name Observer interface
        //@{
        void update() {
            Real v = q_->value();
            QL_REQUIRE(v>0.0, "Invalid quote, cannot take log of non-postive number");
            logValue_= std::log(v);
        }
        //@}
     protected:
        Handle<Quote> q_;
        Real logValue_;
    };

}

#endif
