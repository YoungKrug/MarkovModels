#pragma once
#include <string>

struct SequenceData
{
    enum Type
    {
        CounterClockWise,
        ClockWise
    };
    std::string Sequence;
    Type sequenceType;
    std::string name;
    bool isORF = false;
};
