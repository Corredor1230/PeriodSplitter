#include<vector>
#include<math.h>
#include<stdexcept>
#include"include/SitranoHeader.h"

class ZPFilter{
public:
    ZPFilter(){};
    ~ZPFilter(){};
    static std::vector<float> filtfilt(const std::vector<float>& input, float sr, float cutoffFreq)
    {
        if (input.empty()) return{};
        if (input.size() < 2) return input;

        float u = std::tan(Sihat::PI * cutoffFreq / sr);
        float common = 1.0 + u;

        float b0 = u / common;
        float b1 = b0;
        float a1 = (u - 1.0) / common;

        std::vector<float> output = input;

        int n = input.size();

        float xPrev = input[0];
        float yPrev = input[0];

        for (int i = 0; i < n; ++i)
        {
            float xCurr = input[i];
            float yCurr = (b0 * xCurr) + (b1 * xPrev) - (a1 * yPrev);

            output[i] = yCurr;

            xPrev = xCurr;
            yPrev = yCurr;
        }

        xPrev = output[n - 1];
        yPrev = output[n - 1];

        for (int i = n - 1; i >= 0; --i)
        {
            float xCurr = output[i];
            float yCurr = (b0 * xCurr) + (b1 * xPrev) - (a1 * yPrev);

            output[i] = yCurr;

            xPrev = xCurr;
            yPrev = yCurr;
        }

        return output;
    }  
};