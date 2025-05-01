#ifndef PRICEABLE_H
#define PRICEABLE_H

#include "Black&scholes.h"

class Priceable
{
    public:
        virtual double price(const BlackScholesModel& model) =0; //pure virtual 
}; 


#endif