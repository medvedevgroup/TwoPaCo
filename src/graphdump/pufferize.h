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
constexpr uint16_t MERGE_LEFT = 1 << 11;
constexpr uint16_t MERGE_RIGHT = 1 << 12;
constexpr uint16_t COMPLEX = 1 << 13;

struct KmerInfo {
    //<..>
    //<1 bit: set if kmer has already been seen>
    //<1 bit: set if a reference ends with the kmer in forward>
    //<1 bit: set if a reference starts with the kmer in forward>
    //<4 bits: characters that precede kmer in forward>
    //<4 bits: characters that succeed kmer in forward>
    uint16_t kinf{0};

    void setPrecedingChar(bool isFw, char c) {
        if (isFw) {
            auto idx = TwoPaCo::DnaChar::MakeUpChar(c);
            kinf |= 1 << (4 + idx);
        } else {
            auto idx = TwoPaCo::DnaChar::MakeUpChar(TwoPaCo::DnaChar::ReverseChar(c));
            kinf |= 1 << (idx);
        }
    }
    void setSucceedingChar(bool isFw, char c) {
        if (isFw) {
            auto idx = TwoPaCo::DnaChar::MakeUpChar(c);
            kinf |= 1 << (idx);
        } else {
            auto idx = TwoPaCo::DnaChar::MakeUpChar(TwoPaCo::DnaChar::ReverseChar(c));
            kinf |= 1 << (4 + idx);
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

    void setCropBoth() {
        kinf |= COMPLEX;
    }

    bool cropBoth() {
        return  kinf & COMPLEX;
    }

    void setCropStart() {
        kinf |= MERGE_LEFT;
    }

    bool cropStart() {
        return kinf & MERGE_LEFT;
    }

    void setCropEnd() {
        kinf |= MERGE_RIGHT;
    }

    bool cropEnd() {
        return kinf & MERGE_RIGHT;
    }

    uint16_t countSucceeding() {
        uint16_t cnt{0};
        for (auto i = 0; i < 4; i++) {
            cnt += (kinf & (1 << i));
        }
        return cnt;
    }

    uint16_t countPreceding() {
        uint16_t cnt{0};
        for (auto i = 4; i < 8; i++) {
            cnt += (kinf & (1 << i));
        }
        return cnt;
    }
    void decideType() {
        if (!kinf) // if kmer doesn't exist
            return;
        auto precedeCnt = countPreceding();
        auto succeedCnt = countSucceeding();
        if (precedeCnt > 1 and succeedCnt > 1) {
            setCropBoth();
        } else if (succeedCnt > 1) {
            if (isStart()) {
                setCropBoth();
            } else {
                setCropStart();
            }
        } else if (precedeCnt > 1) {
            if (isEnd()) {
                setCropBoth();
            } else {
                setCropEnd();
            }
        } else if (precedeCnt == 1 and succeedCnt == 1) {
            if (!isEnd() or !isStart()) {
                std::cerr << "ERROR!! Such case cannot happen in the output of TwoPaCo: precedingCnt=1, succeedingCnt=1, neither start nor end\n";
            }
            setCropBoth();
        } // otherwise, we don't require to crop any nucleotides from any sides of a contig/segment

    }
};

#endif //TWOPACO_PUFFERIZE_H
