//
// Created by Fatemeh Almodaresi on 8/3/18.
//

#ifndef TWOPACO_PUFFERIZE_H
#define TWOPACO_PUFFERIZE_H

#include <dnachar.h>
#include <junctionapi/junctionapi.h>

constexpr uint16_t START = 1 << 8;
constexpr uint16_t END = 1 << 9;
constexpr uint16_t SEEN = 1 << 10;

struct KmerInfo {
    //<..>
    //<1 bit: set if kmer has already been seen>
    //<1 bit: set if a reference ends with the kmer in forward>
    //<1 bit: set if a reference starts with the kmer in forward>
    //<4 bits: characters that precede kmer in forward>
    //<4 bits: characters that succeed kmer in forward>
    uint16_t kinf{0};

    void setChar(bool precedes, bool isRC, char c) {
        char trueC = c;
        if (isRC)
            trueC = TwoPaCo::DnaChar::ReverseChar(c);
        auto idx = TwoPaCo::DnaChar::MakeUpChar(c);
        if (precedes == isRC) {
            kinf |= 1 << (idx);
        } else {
            kinf |= 1 << (4+idx);
        }
    }
    void setStart() {
        kinf |= START;
    }
    void setEnd() {
        kinf |= END;
    }
    void setSeen() {
        kinf |= SEEN;
    }

    bool isStart() {
        return kinf & START;
    }

    bool isEnd() {
        return kinf & END;
    }

    bool seen() {
        return kinf & SEEN;
    }

};

#endif //TWOPACO_PUFFERIZE_H
