#define _CRT_SECURE_NO_WARNINGS
 
#include <array>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <cmath>
#include <vector>

#include <algorithm> // I added to fix incompatability with C++11
 
#define _USE_MATH_DEFINES
#include <math.h>
 
//=====================================================================================
// SNumeric - uses phantom types to enforce type safety
//=====================================================================================
template <typename T, typename PHANTOM_TYPE>
struct SNumeric
{
public:
    explicit SNumeric(const T &value) : m_value(value) { }
    SNumeric() : m_value() { }
    inline T& Value() { return m_value; }
    inline const T& Value() const { return m_value; }
 
    typedef SNumeric<T, PHANTOM_TYPE> TType;
    typedef T TInnerType;
 
    // Math Operations
    TType operator+ (const TType &b) const
    {
        return TType(this->Value() + b.Value());
    }
 
    TType operator- (const TType &b) const
    {
        return TType(this->Value() - b.Value());
    }
 
    TType operator* (const TType &b) const
    {
        return TType(this->Value() * b.Value());
    }
 
    TType operator/ (const TType &b) const
    {
        return TType(this->Value() / b.Value());
    }
 
    TType& operator+= (const TType &b)
    {
        Value() += b.Value();
        return *this;
    }
 
    TType& operator-= (const TType &b)
    {
        Value() -= b.Value();
        return *this;
    }
 
    TType& operator*= (const TType &b)
    {
        Value() *= b.Value();
        return *this;
    }
 
    TType& operator/= (const TType &b)
    {
        Value() /= b.Value();
        return *this;
    }
 
    TType& operator++ ()
    {
        Value()++;
        return *this;
    }
 
    TType& operator-- ()
    {
        Value()--;
        return *this;
    }
 
    /* original code:
    // Extended Math Operations
    template <typename T>
    T Divide(const TType &b)
    {
        return ((T)this->Value()) / ((T)b.Value());
    }
    */
    // Extended Math Operations
    template <typename A>
    A Divide(const TType &b)
    {
        return ((A)this->Value()) / ((A)b.Value());
    }
 
    // Logic Operations
    bool operator< (const TType &b) const {
        return this->Value() < b.Value();
    }
    bool operator<= (const TType &b) const {
        return this->Value() <= b.Value();
    }
    bool operator> (const TType &b) const {
        return this->Value() > b.Value();
    }
    bool operator>= (const TType &b) const {
        return this->Value() >= b.Value();
    }
    bool operator== (const TType &b) const {
        return this->Value() == b.Value();
    }
    bool operator!= (const TType &b) const {
        return this->Value() != b.Value();
    }
 
private:
    T m_value;
};
 
//=====================================================================================
// Typedefs
//=====================================================================================
 
typedef uint8_t uint8;
typedef uint16_t uint16;
typedef uint32_t uint32;
typedef int16_t int16;
typedef int32_t int32;
 
// type safe types!
typedef SNumeric<float, struct S__Frequency>      TFrequency;
typedef SNumeric<uint32, struct S__TimeMs>        TTimeMs;
typedef SNumeric<uint32, struct S__Samples>       TSamples;
typedef SNumeric<float, struct S__FractSamples>   TFractionalSamples;
typedef SNumeric<float, struct S__Decibels>       TDecibels;
typedef SNumeric<float, struct S__Amplitude>      TAmplitude;
typedef SNumeric<float, struct S__Phase>          TPhase;
 
//=====================================================================================
// Constants
//=====================================================================================
 
static const float c_pi = (float)M_PI;
static const float c_twoPi = c_pi * 2.0f;
 
//=====================================================================================
// Structs
//=====================================================================================
 
struct SSoundSettings
{
    TSamples        m_sampleRate;
    TTimeMs         m_lengthMs;
    TSamples        m_currentSample;
};
 
//=====================================================================================
// CMultiTapReverb -> the multi tap reverb object
//=====================================================================================
 
struct SReverbTap
{
    TSamples    m_timeOffset;
    TAmplitude  m_feedback;
};
 
class CMultitapReverb
{
public:
    CMultitapReverb(const std::vector<SReverbTap>& taps)
        : m_sampleIndex(0)
    {
        // copy the taps table
        m_taps = taps;
 
        // find out the largest tap time offset so we know how big to make the buffer
        TSamples largestTimeOffset(0);
        std::for_each(m_taps.begin(), m_taps.end(),
            [&largestTimeOffset](const SReverbTap& tap)
            {
                if (tap.m_timeOffset > largestTimeOffset)
                    largestTimeOffset = tap.m_timeOffset;
            }
        );
 
        // if it's 0, bail out, we are done
        if (largestTimeOffset.Value() == 0)
            return;
 
        // else resize our internal buffer and fill it with silence
        m_samples.resize(largestTimeOffset.Value()+1);
        std::fill(m_samples.begin(), m_samples.end(), TAmplitude(0.0f));
    }
 
    TAmplitude ProcessSample(TAmplitude sample)
    {
        // if no taps, or none with any time value, bail out!
        if (m_samples.size() == 0)
            return sample;
 
        // take our taps from the delay buffer
        TAmplitude outSample = sample;
        std::for_each(m_taps.begin(), m_taps.end(),
            [&outSample, this](const SReverbTap& tap)
            {
                size_t tapSampleIndex;
                if (tap.m_timeOffset.Value() > m_sampleIndex)
                    tapSampleIndex = m_samples.size() - 1 - (tap.m_timeOffset.Value() - m_sampleIndex);
                else
                    tapSampleIndex = m_sampleIndex - tap.m_timeOffset.Value();
 
                outSample += m_samples[tapSampleIndex] * tap.m_feedback;
            }
        );
 
        // put the output sample into the buffer
        m_samples[m_sampleIndex] = outSample;
 
        // advance the circular buffer index
        m_sampleIndex++;
        if (m_sampleIndex >= m_samples.size())
            m_sampleIndex = 0;
 
        // return the reverbed sample
        return outSample;
    }
 
private:
    std::vector<SReverbTap> m_taps;
    std::vector<TAmplitude> m_samples;
    size_t                  m_sampleIndex;
};
 
//=====================================================================================
// Conversion Functions
//=====================================================================================
inline TDecibels AmplitudeToDB(TAmplitude volume)
{
    return TDecibels(log10(volume.Value()));
}
 
inline TAmplitude DBToAmplitude(TDecibels dB)
{
    return TAmplitude(pow(10.0f, dB.Value() / 20.0f));
}
 
TSamples SecondsToSamples(const SSoundSettings &s, float seconds)
{
    return TSamples((int)(seconds * (float)s.m_sampleRate.Value()));
}
 
TSamples MilliSecondsToSamples(const SSoundSettings &s, float milliseconds)
{
    return SecondsToSamples(s, milliseconds / 1000.0f);
}
 
TTimeMs SecondsToMilliseconds(float seconds)
{
    return TTimeMs((uint32)(seconds * 1000.0f));
}
 
TFrequency Frequency(float octave, float note)
{
    /* frequency = 440×(2^(n/12))
    Notes:
    0  = A
    1  = A#
    2  = B
    3  = C
    4  = C#
    5  = D
    6  = D#
    7  = E
    8  = F
    9  = F#
    10 = G
    11 = G# */
    return TFrequency((float)(440 * pow(2.0, ((double)((octave - 4) * 12 + note)) / 12.0)));
}
 
template <typename T>
T AmplitudeToAudioSample(const TAmplitude& in)
{
    const T c_min = std::numeric_limits<T>::min();
    const T c_max = std::numeric_limits<T>::max();
    const float c_minFloat = (float)c_min;
    const float c_maxFloat = (float)c_max;
 
    float ret = in.Value() * c_maxFloat;
 
    if (ret < c_minFloat)
        return c_min;
 
    if (ret > c_maxFloat)
        return c_max;
 
    return (T)ret;
}
 
TAmplitude GetLerpedAudioSample(const std::vector<TAmplitude>& samples, TFractionalSamples& index)
{
    // get the index of each sample and the fractional blend amount
    uint32 a = (uint32)floor(index.Value());
    uint32 b = a + 1;
    float fract = index.Value() - floor(index.Value());
 
    // get our two amplitudes
    float ampA = 0.0f;
    if (a >= 0 && a < samples.size())
        ampA = samples[a].Value();
 
    float ampB = 0.0f;
    if (b >= 0 && b < samples.size())
        ampB = samples[b].Value();
 
    // return the lerped result
    return TAmplitude(fract * ampB + (1.0f - fract) * ampA);
}
 
void NormalizeSamples(std::vector<TAmplitude>& samples, TAmplitude maxAmplitude)
{
    // nothing to do if no samples
    if (samples.size() == 0)
        return;
 
    // 1) find the largest absolute value in the samples.
    TAmplitude largestAbsVal = TAmplitude(abs(samples.front().Value()));
    std::for_each(samples.begin() + 1, samples.end(), [&largestAbsVal](const TAmplitude &a)
    {
        TAmplitude absVal = TAmplitude(abs(a.Value()));
        if (absVal > largestAbsVal)
            largestAbsVal = absVal;
    }
    );
 
    // 2) adjust largestAbsVal so that when we divide all samples, none will be bigger than maxAmplitude
    // if the value we are going to divide by is <= 0, bail out
    largestAbsVal /= maxAmplitude;
    if (largestAbsVal <= TAmplitude(0.0f))
        return;
 
    // 3) divide all numbers by the largest absolute value seen so all samples are [-maxAmplitude,+maxAmplitude]
    std::for_each(samples.begin(), samples.end(), [&largestAbsVal](TAmplitude &a)
    {
        a /= largestAbsVal;
 
        if (a >= TAmplitude(1.0f))
        {
            int ijkl = 0;
        }
    }
    );
}
 
void ResampleData(std::vector<TAmplitude>& samples, int srcSampleRate, int destSampleRate)
{
    //if the requested sample rate is the sample rate it already is, bail out and do nothing
    if (srcSampleRate == destSampleRate)
        return;
 
    //calculate the ratio of the old sample rate to the new
    float fResampleRatio = (float)destSampleRate / (float)srcSampleRate;
 
    //calculate how many samples the new data will have and allocate the new sample data
    int nNewDataNumSamples = (int)((float)samples.size() * fResampleRatio);
 
    std::vector<TAmplitude> newSamples;
    newSamples.resize(nNewDataNumSamples);
 
    //get each lerped output sample.  There are higher quality ways to resample
    for (int nIndex = 0; nIndex < nNewDataNumSamples; ++nIndex) {
        auto frac = TFractionalSamples((float)nIndex / fResampleRatio);
        newSamples[nIndex] = GetLerpedAudioSample(samples, frac);
    }
 
    //free the old data and set the new data
    std::swap(samples, newSamples);
}
 
void ChangeNumChannels(std::vector<TAmplitude>& samples, int nSrcChannels, int nDestChannels)
{
    //if the number of channels requested is the number of channels already there, or either number of channels is not mono or stereo, return
    if (nSrcChannels == nDestChannels ||
        nSrcChannels < 1 || nSrcChannels > 2 ||
        nDestChannels < 1 || nDestChannels > 2)
    {
        return;
    }
 
    //if converting from mono to stereo, duplicate the mono channel to make stereo
    if (nDestChannels == 2)
    {
        std::vector<TAmplitude> newSamples;
        newSamples.resize(samples.size() * 2);
        for (size_t index = 0; index < samples.size(); ++index)
        {
            newSamples[index * 2] = samples[index];
            newSamples[index * 2 + 1] = samples[index];
        }
 
        std::swap(samples, newSamples);
    }
    //else converting from stereo to mono, mix the stereo channels together to make mono
    else
    {
        std::vector<TAmplitude> newSamples;
        newSamples.resize(samples.size() / 2);
        for (size_t index = 0; index < samples.size() / 2; ++index)
            newSamples[index] = samples[index * 2] + samples[index * 2 + 1];
 
        std::swap(samples, newSamples);
    }
}
 
float PCMToFloat(unsigned char *pPCMData, int nNumBytes)
{
    switch (nNumBytes)
    {
    case 1:
    {
        uint8 data = pPCMData[0];
        return (float)data / 255.0f;
    }
    case 2:
    {
        int16 data = pPCMData[1] << 8 | pPCMData[0];
        return ((float)data) / ((float)0x00007fff);
    }
    case 3:
    {
        int32 data = pPCMData[2] << 16 | pPCMData[1] << 8 | pPCMData[0];
        return ((float)data) / ((float)0x007fffff);
    }
    case 4:
    {
        int32 data = pPCMData[3] << 24 | pPCMData[2] << 16 | pPCMData[1] << 8 | pPCMData[0];
        return ((float)data) / ((float)0x7fffffff);
    }
    default:
    {
        return 0.0f;
    }
    }
}
 
//=====================================================================================
// Wave File Writing Code
//=====================================================================================
struct SMinimalWaveFileHeader
{
    //the main chunk
    unsigned char m_szChunkID[4];      //0
    uint32        m_nChunkSize;        //4
    unsigned char m_szFormat[4];       //8
 
    //sub chunk 1 "fmt "
    unsigned char m_szSubChunk1ID[4];  //12
    uint32        m_nSubChunk1Size;    //16
    uint16        m_nAudioFormat;      //18
    uint16        m_nNumChannels;      //20
    uint32        m_nSampleRate;       //24
    uint32        m_nByteRate;         //28
    uint16        m_nBlockAlign;       //30
    uint16        m_nBitsPerSample;    //32
 
    //sub chunk 2 "data"
    unsigned char m_szSubChunk2ID[4];  //36
    uint32        m_nSubChunk2Size;    //40
 
    //then comes the data!
};
 
//this writes a wave file
template <typename T>
bool WriteWaveFile(const char *fileName, const std::vector<TAmplitude> &samples, const SSoundSettings &sound)
{
    //open the file if we can
    FILE *file = fopen(fileName, "w+b");
    if (!file)
        return false;
 
    //calculate bits per sample and the data size
    const int32 bitsPerSample = sizeof(T) * 8;
    const int dataSize = samples.size() * sizeof(T);
 
    SMinimalWaveFileHeader waveHeader;
 
    //fill out the main chunk
    memcpy(waveHeader.m_szChunkID, "RIFF", 4);
    waveHeader.m_nChunkSize = dataSize + 36;
    memcpy(waveHeader.m_szFormat, "WAVE", 4);
 
    //fill out sub chunk 1 "fmt "
    memcpy(waveHeader.m_szSubChunk1ID, "fmt ", 4);
    waveHeader.m_nSubChunk1Size = 16;
    waveHeader.m_nAudioFormat = 1;
    waveHeader.m_nNumChannels = 1;
    waveHeader.m_nSampleRate = sound.m_sampleRate.Value();
    waveHeader.m_nByteRate = sound.m_sampleRate.Value() * 1 * bitsPerSample / 8;
    waveHeader.m_nBlockAlign = 1 * bitsPerSample / 8;
    waveHeader.m_nBitsPerSample = bitsPerSample;
 
    //fill out sub chunk 2 "data"
    memcpy(waveHeader.m_szSubChunk2ID, "data", 4);
    waveHeader.m_nSubChunk2Size = dataSize;
 
    //write the header
    fwrite(&waveHeader, sizeof(SMinimalWaveFileHeader), 1, file);
 
    //write the wave data itself, converting it from float to the type specified
    std::vector<T> outSamples;
    outSamples.resize(samples.size());
    for (size_t index = 0; index < samples.size(); ++index)
        outSamples[index] = AmplitudeToAudioSample<T>(samples[index]);
    fwrite(&outSamples[0], dataSize, 1, file);
 
    //close the file and return success
    fclose(file);
    return true;
}
 
//loads a wave file in.  Converts from source format into the specified format
// TOTAL HONESTY: some wave files seem to have problems being loaded through this function and I don't have
// time to investigate why.  It seems to work best with 16 bit mono wave files.
// If you need more robust file loading, check out libsndfile at http://www.mega-nerd.com/libsndfile/
bool ReadWaveFile(const char *fileName, std::vector<TAmplitude>& samples, int32 sampleRate)
{
    //open the file if we can
    FILE *File = fopen(fileName, "rb");
    if (!File)
    {
        return false;
    }
 
    //read the main chunk ID and make sure it's "RIFF"
    char buffer[5];
    buffer[4] = 0;
    if (fread(buffer, 4, 1, File) != 1 || strcmp(buffer, "RIFF"))
    {
        fclose(File);
        return false;
    }
 
    //read the main chunk size
    uint32 nChunkSize;
    if (fread(&nChunkSize, 4, 1, File) != 1)
    {
        fclose(File);
        return false;
    }
 
    //read the format and make sure it's "WAVE"
    if (fread(buffer, 4, 1, File) != 1 || strcmp(buffer, "WAVE"))
    {
        fclose(File);
        return false;
    }
 
    long chunkPosFmt = -1;
    long chunkPosData = -1;
 
    while (chunkPosFmt == -1 || chunkPosData == -1)
    {
        //read a sub chunk id and a chunk size if we can
        if (fread(buffer, 4, 1, File) != 1 || fread(&nChunkSize, 4, 1, File) != 1)
        {
            fclose(File);
            return false;
        }
 
        //if we hit a fmt
        if (!strcmp(buffer, "fmt "))
        {
            chunkPosFmt = ftell(File) - 8;
        }
        //else if we hit a data
        else if (!strcmp(buffer, "data"))
        {
            chunkPosData = ftell(File) - 8;
        }
 
        //skip to the next chunk
        fseek(File, nChunkSize, SEEK_CUR);
    }
 
    //we'll use this handy struct to load in 
    SMinimalWaveFileHeader waveData;
 
    //load the fmt part if we can
    fseek(File, chunkPosFmt, SEEK_SET);
    if (fread(&waveData.m_szSubChunk1ID, 24, 1, File) != 1)
    {
        fclose(File);
        return false;
    }
 
    //load the data part if we can
    fseek(File, chunkPosData, SEEK_SET);
    if (fread(&waveData.m_szSubChunk2ID, 8, 1, File) != 1)
    {
        fclose(File);
        return false;
    }
 
    //verify a couple things about the file data
    if (waveData.m_nAudioFormat != 1 ||       //only pcm data
        waveData.m_nNumChannels < 1 ||        //must have a channel
        waveData.m_nNumChannels > 2 ||        //must not have more than 2
        waveData.m_nBitsPerSample > 32 ||     //32 bits per sample max
        waveData.m_nBitsPerSample % 8 != 0 || //must be a multiple of 8 bites
        waveData.m_nBlockAlign > 8)           //blocks must be 8 bytes or lower
    {
        fclose(File);
        return false;
    }
 
    //figure out how many samples and blocks there are total in the source data
    int nBytesPerBlock = waveData.m_nBlockAlign;
    int nNumBlocks = waveData.m_nSubChunk2Size / nBytesPerBlock;
    int nNumSourceSamples = nNumBlocks * waveData.m_nNumChannels;
 
    //allocate space for the source samples
    samples.resize(nNumSourceSamples);
 
    //maximum size of a block is 8 bytes.  4 bytes per samples, 2 channels
    unsigned char pBlockData[8];
    memset(pBlockData, 0, 8);
 
    //read in the source samples at whatever sample rate / number of channels it might be in
    int nBytesPerSample = nBytesPerBlock / waveData.m_nNumChannels;
    for (int nIndex = 0; nIndex < nNumSourceSamples; nIndex += waveData.m_nNumChannels)
    {
        //read in a block
        if (fread(pBlockData, waveData.m_nBlockAlign, 1, File) != 1)
        {
            fclose(File);
            return false;
        }
 
        //get the first sample
        samples[nIndex].Value() = PCMToFloat(pBlockData, nBytesPerSample);
 
        //get the second sample if there is one
        if (waveData.m_nNumChannels == 2)
        {
            samples[nIndex + 1].Value() = PCMToFloat(&pBlockData[nBytesPerSample], nBytesPerSample);
        }
    }
 
    //re-sample the sample rate up or down as needed
    ResampleData(samples, waveData.m_nSampleRate, sampleRate);
 
    //handle switching from mono to stereo or vice versa
    ChangeNumChannels(samples, waveData.m_nNumChannels, 1);
 
    return true;
}
 
//=====================================================================================
// Oscilators
//=====================================================================================
 
void AdvancePhase(TPhase &phase, TFrequency frequency, TSamples sampleRate)
{
    phase += TPhase(frequency.Value() / (float)sampleRate.Value());
    while (phase >= TPhase(1.0f))
        phase -= TPhase(1.0f);
    while (phase < TPhase(0.0f))
        phase += TPhase(1.0f);
}
 
TAmplitude AdvanceOscilator_Sine(TPhase &phase, TFrequency frequency, TSamples sampleRate)
{
    AdvancePhase(phase, frequency, sampleRate);
    return TAmplitude(sin(phase.Value()*c_twoPi));
}
 
TAmplitude AdvanceOscilator_Saw(TPhase &phase, TFrequency frequency, TSamples sampleRate)
{
    AdvancePhase(phase, frequency, sampleRate);
    return TAmplitude(phase.Value() * 2.0f - 1.0f);
}
 
TAmplitude AdvanceOscilator_Square(TPhase &phase, TFrequency frequency, TSamples sampleRate)
{
    AdvancePhase(phase, frequency, sampleRate);
    return TAmplitude(phase.Value() > 0.5f ? 1.0f : -1.0f);
}
 
TAmplitude AdvanceOscilator_Triangle(TPhase &phase, TFrequency frequency, TSamples sampleRate)
{
    AdvancePhase(phase, frequency, sampleRate);
    if (phase > TPhase(0.5f))
        return TAmplitude((((1.0f - phase.Value()) * 2.0f) * 2.0f) - 1.0f);
    else
        return TAmplitude(((phase.Value() * 2.0f) * 2.0f) - 1.0f);
}
 
TAmplitude AdvanceOscilator_Saw_BandLimited(TPhase &phase, TFrequency frequency, TSamples sampleRate)
{
    AdvancePhase(phase, frequency, sampleRate);
 
    // sum the harmonics
    TAmplitude ret(0.0f);
    for (int harmonicIndex = 1; harmonicIndex <= 4; ++harmonicIndex)
    {
        TPhase harmonicPhase = phase * TPhase((float)harmonicIndex);
        ret += TAmplitude(sin(harmonicPhase.Value()*c_twoPi) / (float)harmonicIndex);
    }
 
    //adjust the volume
    ret *= TAmplitude(2.0f / c_pi);
 
    return ret;
}
 
TAmplitude AdvanceOscilator_Square_BandLimited(TPhase &phase, TFrequency frequency, TSamples sampleRate)
{
    AdvancePhase(phase, frequency, sampleRate);
 
    // sum the harmonics
    TAmplitude ret(0.0f);
    for (int harmonicIndex = 1; harmonicIndex <= 4; ++harmonicIndex)
    {
        float harmonicFactor = (float)harmonicIndex * 2.0f - 1.0f;
        TPhase harmonicPhase = phase * TPhase(harmonicFactor);
        ret += TAmplitude(sin(harmonicPhase.Value()*c_twoPi) / harmonicFactor);
    }
 
    //adjust the volume
    ret *= TAmplitude(4.0f / c_pi);
 
    return ret;
}
 
TAmplitude AdvanceOscilator_Triangle_BandLimited(TPhase &phase, TFrequency frequency, TSamples sampleRate)
{
    AdvancePhase(phase, frequency, sampleRate);
 
    // sum the harmonics
    TAmplitude ret(0.0f);
    TAmplitude signFlip(1.0f);
    for (int harmonicIndex = 1; harmonicIndex <= 4; ++harmonicIndex)
    {
        float harmonicFactor = (float)harmonicIndex * 2.0f - 1.0f;
        TPhase harmonicPhase = phase * TPhase(harmonicFactor);
        ret += TAmplitude(sin(harmonicPhase.Value()*c_twoPi) / (harmonicFactor*harmonicFactor)) * signFlip;
        signFlip *= TAmplitude(-1.0f);
    }
 
    //adjust the volume
    ret *= TAmplitude(8.0f / (c_pi*c_pi));
 
    return ret;
}
 
//=====================================================================================
// Main
//=====================================================================================
int main(int argc, char **argv)
{
    //our desired sound parameters
    SSoundSettings sound;
    sound.m_sampleRate = TSamples(44100);
    sound.m_lengthMs = SecondsToMilliseconds(4.0f);
 
    // create a reverb object with a list of taps
    CMultitapReverb reverb(
        {
            { MilliSecondsToSamples(sound,  79.0f), DBToAmplitude(TDecibels(-25.0f)) },
            { MilliSecondsToSamples(sound, 130.0f), DBToAmplitude(TDecibels(-23.0f)) },
            { MilliSecondsToSamples(sound, 230.0f), DBToAmplitude(TDecibels(-15.0f)) },
            { MilliSecondsToSamples(sound, 340.0f), DBToAmplitude(TDecibels(-23.0f)) },
            { MilliSecondsToSamples(sound, 470.0f), DBToAmplitude(TDecibels(-17.0f)) },
            { MilliSecondsToSamples(sound, 532.0f), DBToAmplitude(TDecibels(-21.0f)) },
            { MilliSecondsToSamples(sound, 662.0f), DBToAmplitude(TDecibels(-13.0f)) },
        }
    );
 
    // load the wave file if we can
    std::vector<TAmplitude> inputData;
    if (!ReadWaveFile("in.wav", inputData, sound.m_sampleRate.Value()))
    {
        printf("could not load wave file!");
        return 0;
    }
 
    // allocate space for the output file
    std::vector<TAmplitude> samples;
    samples.resize(inputData.size());
 
    //apply the delay effect to the file
    const TSamples c_envelopeSize = MilliSecondsToSamples(sound, 50.0f);
    for (TSamples index = TSamples(0), numSamples(samples.size()); index < numSamples; ++index)
    {
        // calculate envelope at front and end of sound.
        TAmplitude envelope(1.0f);
        if (index < c_envelopeSize)
            envelope = TAmplitude((float)index.Value() / (float)c_envelopeSize.Value());
        else if (index >(numSamples - c_envelopeSize))
            envelope = TAmplitude(1.0f) - TAmplitude((float)(index - (numSamples - c_envelopeSize)).Value() / (float)c_envelopeSize.Value());
 
        // put our input through the reverb process
        TAmplitude outSample = reverb.ProcessSample(inputData[index.Value()]);
 
        // mix the sample with the offset sample and apply the envelope for the front and back of the sound
        samples[index.Value()] = outSample * envelope;
    }
 
    // normalize the amplitude of the samples to make sure they are as loud as possible without clipping
    // give 3db of headroom
    NormalizeSamples(samples, DBToAmplitude(TDecibels(-3.0f)));
 
    // save as a wave file
    WriteWaveFile<int16_t>("out.wav", samples, sound);
 
    return 0;
}
