#pragma once

#include<vector>
#include<algorithm>
#include<cmath>
#include<string>
namespace SihatInterpolation
{
    enum class Type : uint32_t
    {
        Linear = 0,
        Spline = 1,
        Hybrid = 2
    };

    inline float f_interpolate(const float valA, const float valB, const float alpha, const Type iType = Type::Linear)
    {
        float outVal = 0.0;

        switch(iType)
        {
            case Type::Linear:
            {
                outVal = valA + alpha * (valB - valA);
                break;
            }
            case Type::Spline:
            {
                break;
            }
            default:
            {
                outVal = valA + alpha * (valB - valA);
                break;
            }
        }

        return outVal;
    }

    inline int i_interpolate(const int valA, const int valB, const float alpha, const Type iType = Type::Linear)
    {
        int outVal = 0;

        switch(iType)
        {
            case Type::Linear:
            {
                float exactVal = valA + alpha * (valB - valA);
                outVal = static_cast<int>(std::round(exactVal));
                break; 
            }
            case Type::Spline:
            {
                break; 
            }
            default:
            {
                float exactVal = valA + alpha * (valB - valA);
                outVal = static_cast<int>(std::round(exactVal));
                break;
            }
        }

        return outVal;
    }

    inline std::vector<float> fvec_interpolate(const std::vector<float>& vecA, const std::vector<float>& vecB, const float alpha, const Type iType)
    {
        if (iType != Type::Linear) 
        {
            return {}; // Or handle your other InterpTypes here
        }

        const size_t sizeA = vecA.size();
        const size_t sizeB = vecB.size();
        const size_t minSize = std::min(sizeA, sizeB);
        const size_t maxSize = std::max(sizeA, sizeB);

        std::vector<float> out;
        out.reserve(maxSize); // Prevent unnecessary reallocations
        float iVal = 0.0;

        for (size_t i = 0; i < minSize; ++i)
        {
            iVal = f_interpolate(vecA[i], vecB[i], alpha, iType);
            out.push_back(iVal);
        }

        if (sizeA > sizeB)
        {
            for (size_t i = minSize; i < maxSize; ++i)
            {
                iVal = f_interpolate(vecA[i], 0.0, alpha, iType);
                out.push_back(iVal);
            }
        }
        else if (sizeB > sizeA)
        {
            for (size_t i = minSize; i < maxSize; ++i)
            {
                iVal = f_interpolate(0.0, vecB[i], alpha, iType);
                out.push_back(iVal); 
            }
        }

        return out;
    }

    inline std::vector<int> ivec_interpolate(const std::vector<int>& vecA, const std::vector<int>& vecB, const float alpha, const Type iType)
    {
        if (iType != Type::Linear) 
        {
            return {}; // Or handle your other InterpTypes here
        }

        const size_t sizeA = vecA.size();
        const size_t sizeB = vecB.size();
        const size_t minSize = std::min(sizeA, sizeB);
        const size_t maxSize = std::max(sizeA, sizeB);

        std::vector<int> out;
        out.reserve(maxSize); // Prevent unnecessary reallocations
        float iVal = 0.0;

        for (size_t i = 0; i < minSize; ++i)
        {
            iVal = i_interpolate(vecA[i], vecB[i], alpha, iType);
            out.push_back(iVal);
        }

        if (sizeA > sizeB)
        {
            for (size_t i = minSize; i < maxSize; ++i)
            {
                iVal = i_interpolate(vecA[i], 0.0, alpha, iType);
                out.push_back(iVal);
            }
        }
        else if (sizeB > sizeA)
        {
            for (size_t i = minSize; i < maxSize; ++i)
            {
                iVal = i_interpolate(0.0, vecB[i], alpha, iType);
                out.push_back(iVal); 
            }
        }

        return out;
    }

    inline std::vector<uint32_t> ivec_interpolate(const std::vector<uint32_t>& vecA, const std::vector<uint32_t>& vecB, const float alpha, const Type iType)
    {
        if (iType != Type::Linear) 
        {
            return {}; // Or handle your other InterpTypes here
        }

        const size_t sizeA = vecA.size();
        const size_t sizeB = vecB.size();
        const size_t minSize = std::min(sizeA, sizeB);
        const size_t maxSize = std::max(sizeA, sizeB);

        std::vector<uint32_t> out;
        out.reserve(maxSize); // Prevent unnecessary reallocations
        float iVal = 0.0;

        for (size_t i = 0; i < minSize; ++i)
        {
            iVal = i_interpolate(vecA[i], vecB[i], alpha, iType);
            out.push_back(iVal);
        }

        if (sizeA > sizeB)
        {
            for (size_t i = minSize; i < maxSize; ++i)
            {
                iVal = i_interpolate(vecA[i], 0.0, alpha, iType);
                out.push_back(iVal);
            }
        }
        else if (sizeB > sizeA)
        {
            for (size_t i = minSize; i < maxSize; ++i)
            {
                iVal = i_interpolate(0.0, vecB[i], alpha, iType);
                out.push_back(iVal); 
            }
        }

        return out;
    }

}