#ifndef RNAP_H_
#define RNAP_H_

#include "utils.h"
#include "Transcript.h"

class RNAP 
{
    public:
        RNAP(void); // default constructor
        RNAP(const Transcript& t);
        uint    tr_id_;
        int     strand_;
        DNApos  pos_;
        DNApos  last_pos_;
        void    move(void);
        bool    hasfinished(void);
};


#endif //RNAP_H_