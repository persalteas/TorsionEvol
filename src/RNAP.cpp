#include "RNAP.h"

RNAP::RNAP(void)
{
    strand_ = 0;
    pos_ = -1;
    last_pos_ = -1;
    tr_id_ = -1;
}

RNAP::RNAP(const Transcript& t)
{
    tr_id_ = t.TUindex_;
    strand_ = t.s_;
    pos_ = t.start_;
    last_pos_ = t.end_;
}

void RNAP::move(void)
{
    pos_ += strand_;
}

bool RNAP::hasfinished(void)
{
    return (pos_==last_pos_);
}