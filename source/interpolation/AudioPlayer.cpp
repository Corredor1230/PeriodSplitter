#include "interpolation/AudioPlayer.h"

AudioPlayer::AudioPlayer() {
    if (dac.getDeviceCount() < 1) {
        std::cerr << "[Audio] no audio devices found! \n";
        return;
    }
}

AudioPlayer::~AudioPlayer() {
    if (dac.isStreamOpen()) {
        dac.closeStream();
    }
}

bool AudioPlayer::init(uint32_t sampleRate) {
    RtAudio::StreamParameters parameters;
    parameters.deviceId = dac.getDefaultOutputDevice();
    parameters.nChannels = 1;
    parameters.firstChannel = 0;

    uint32_t bufferFrames = 512;

    try {
        dac.openStream(&parameters, nullptr, RTAUDIO_FLOAT32, sampleRate, &bufferFrames, &AudioPlayer::audioCallback, this);
        dac.startStream();
        return true;
    }
    catch (const std::exception& e) {
        std::cerr << "[Audio Error] " << e.what() << '\n';
        return false;
    }
}

void AudioPlayer::playAudio(const std::vector<float>& newAudio) {
    isPlaying = false;
    audioData = newAudio;
    playhead = 0;
    isPlaying = true;
}

void AudioPlayer::stopAudio() {
    isPlaying = false;
}