#pragma once

#include<RtAudio.h>
#include<vector>
#include<atomic>
#include<iostream>

class AudioPlayer {
public:
    AudioPlayer();
    ~AudioPlayer();

    bool init(uint32_t sampleRate = 96000);
    void playAudio(const std::vector<float>& newAudio);
    void stopAudio();

private:
    RtAudio dac;
    std::vector<float> audioData;
    std::atomic<size_t> playhead{0};
    std::atomic<bool> isPlaying{false};

    static int audioCallback(void* outputBuffer, void* /*inputBuffer*/, unsigned int nBufferFrames,
                             double /*streamTime*/, RtAudioStreamStatus /*status*/, void* userData)
    {
        AudioPlayer* player = static_cast<AudioPlayer*>(userData);
        float* out = static_cast<float*>(outputBuffer);

        // If stopped, just output silence
        if (!player->isPlaying) {
            std::fill(out, out + nBufferFrames, 0.0f);
            return 0;
        }

        size_t current = player->playhead.load();
        size_t total = player->audioData.size();

        for (unsigned int i = 0; i < nBufferFrames; i++) {
            if (current < total) {
                out[i] = player->audioData[current++];
            } else {
                out[i] = 0.0f; // End of file
                player->isPlaying = false;
            }
        }

        player->playhead.store(current);
        return 0;
    }


};