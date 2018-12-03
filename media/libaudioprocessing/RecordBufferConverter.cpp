/*
 * Copyright (C) 2017 The Android Open Source Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#define LOG_TAG "RecordBufferConverter"
//#define LOG_NDEBUG 0

#define FILTER_PI 3.141592653589793

#include <audio_utils/primitives.h>
#include <audio_utils/format.h>
#include <media/AudioMixer.h>  // for UNITY_GAIN_FLOAT
#include <media/AudioResampler.h>
#include <media/BufferProviders.h>
#include <media/RecordBufferConverter.h>
#include <utils/Log.h>

#ifndef ARRAY_SIZE
#define ARRAY_SIZE(x) (sizeof(x)/sizeof((x)[0]))
#endif

template <typename T>
static inline T max(const T& a, const T& b)
{
    return a > b ? a : b;
}

namespace android {

RecordBufferConverter::RecordBufferConverter(
        audio_channel_mask_t srcChannelMask, audio_format_t srcFormat,
        uint32_t srcSampleRate,
        audio_channel_mask_t dstChannelMask, audio_format_t dstFormat,
        uint32_t dstSampleRate) :
            mSrcChannelMask(AUDIO_CHANNEL_INVALID), // updateParameters will set following vars
            // mSrcFormat
            // mSrcSampleRate
            // mDstChannelMask
            // mDstFormat
            // mDstSampleRate
            // mSrcChannelCount
            // mDstChannelCount
            // mDstFrameSize
            mBuf(NULL), mBufFrames(0), mBufFrameSize(0),
            mResampler(NULL),
            mIsLegacyDownmix(false),
            mIsLegacyUpmix(false),
            mRequiresFloat(false),
            mInputConverterProvider(NULL)
{
    (void)updateParameters(srcChannelMask, srcFormat, srcSampleRate,
            dstChannelMask, dstFormat, dstSampleRate);
}

RecordBufferConverter::~RecordBufferConverter() {
    free(mBuf);
    delete mResampler;
    delete mInputConverterProvider;
}

void RecordBufferConverter::apply_low_pass(void *buffer, size_t frames) {
    //ALOGE("TESTFILTER start of apply_low_pass");
    //ALOGE("TESTFILTER mSrcSampleRate: %u, mDstSampleRate: %u", mSrcSampleRate, mDstSampleRate);
    //ALOGE("TESTFILTER mSrcFormat: %#x, mDstFormat: %#x", mSrcFormat, mDstFormat);
    //ALOGE("TESTFILTER mBufFrames: %u, mBufFrameSize: %#x", mBufFrames, mBufFrameSize);
    //ALOGE("TESTFILTER mSrcChannelCount: %u, mDstChannelCount: %u", mSrcChannelCount, mDstChannelCount);
    //ALOGE("TESTFILTER audio_bytes_per_sample(mSrcFormat): %u", audio_bytes_per_sample(mSrcFormat));
    //ALOGE("TESTFILTER audio_bytes_per_sample(mDstFormat): %u", audio_bytes_per_sample(mDstFormat));
    //ALOGE("TESTFILTER frames: %u", frames);

    if (!mFilterIsInited) {
        ALOGE("TESTFILTER initing filter.");

        // Fs: sample rate
        // F0: the cutting frequency
        // Q: attenuation in Db

        uint32_t fs;
        double f0 = 0;
        double dbGain = 0;
        double Q = 0;
        double w0 = 0;
        double c = 0;
        double s = 0;
        double alpha = 0;
        double A = 0;
        double As = 0;
        double a0 = 0;

        if (0) { // low pass
            // args for filter
            fs = mDstSampleRate;
            f0 = 17000;
            Q = 80.0;

            // precomputation
            w0 = 2*FILTER_PI*f0/fs;
            c = cos(w0);
            s = sin(w0);
            alpha = s/(2*Q);

            // parameters
            a0 = 1 + alpha;
            mFilterA1 = -2 * c; mFilterA1 /= a0;
            mFilterA2 = 1 - alpha; mFilterA2 /= a0;
            mFilterB0 = (1 - c)/2; mFilterB0 /= a0;
            mFilterB1 = 1 - c; mFilterB1 /= a0;
            mFilterB2 = (1 - c)/2; mFilterB2 /= a0;
        } else if (1) { // high shelf
            // args for filter
            fs = mDstSampleRate;
            f0 = 17500; // 18000
            Q = 1; // S = 1
            dbGain = -120.0; // -70

            // precomputation
            A = pow(10, dbGain/40);
            As = sqrt(A);
            w0 = 2*FILTER_PI*f0/fs;
            c = cos(w0);
            s = sin(w0);
            alpha = s/(2*Q);

            // parameters
            a0 = (A+1) - (A-1)*c + 2*As*alpha;
            mFilterA1 = 2*( (A-1) - (A+1)*c ); mFilterA1 /= a0;
            mFilterA2 = (A+1) - (A-1)*c - 2*As*alpha; mFilterA2 /= a0;
            mFilterB0 = A*( (A+1) + (A-1)*c + 2*As*alpha ); mFilterB0 /= a0;
            mFilterB1 = -2*A*( (A-1) + (A+1)*c ); mFilterB1 /= a0;
            mFilterB2 = A*( (A+1) + (A-1)*c -2*As*alpha ); mFilterB2 /= a0;
        } else {
            ALOGE("TESTFILTER ERROR no filter specified");
        }

        ALOGE("TESTFILTER fs: %u", fs);
        ALOGE("TESTFILTER f0: %f", f0);
        ALOGE("TESTFILTER dbGain: %G", dbGain);
        ALOGE("TESTFILTER Q: %G", Q);
        ALOGE("TESTFILTER w0: %G", w0);
        ALOGE("TESTFILTER c: %G", c);
        ALOGE("TESTFILTER s: %G", s);
        ALOGE("TESTFILTER alpha: %G", alpha);
        ALOGE("TESTFILTER A: %G", A);
        ALOGE("TESTFILTER As: %G", As);
        ALOGE("TESTFILTER a0: %G", a0);
        ALOGE("TESTFILTER mFilterB0: %G", mFilterB0);
        ALOGE("TESTFILTER mFilterB1: %G", mFilterB1);
        ALOGE("TESTFILTER mFilterB2: %G", mFilterB2);
        ALOGE("TESTFILTER mFilterA1: %G", mFilterA1);
        ALOGE("TESTFILTER mFilterA2: %G", mFilterA2);

        mFilterX1 = 0;
        mFilterX2 = 0;
        mFilterY1 = 0;
        mFilterY2 = 0;

        mFilterIsInited = true;
    }

    size_t frames_num = frames * mDstChannelCount;
    signed short *buf = (signed short *) buffer; // FIXME this should depend re. format size
    //unsigned char *buf_bytes = (unsigned char *) buffer;
    for (size_t i = 0; i < frames_num; i++) {
        // for debugging
        //signed short in_us = buf[i];
        //unsigned char in_byte_0 = buf_bytes[i];
        //unsigned char in_byte_1 = buf_bytes[i+1];

        double x = (double) buf[i];
        double y = mFilterB0 * x + mFilterB1 * mFilterX1 + mFilterB2 * mFilterX2 - mFilterA1 * mFilterY1 - mFilterA2 * mFilterY2;
        mFilterX2 = mFilterX1;
        mFilterX1 = x;
        mFilterY2 = mFilterY1;
        mFilterY1 = y;
        buf[i] = (signed short) (y*0.90);

        //buf[i] = sin(((double) (1000.0 * (2 * FILTER_PI) * i)) / mDstSampleRate) * 30000; // signal at 1KHz
        //buf[i] = (signed short) (((double) buf[i]) * 2);
        //ALOGE("TESTFILTER values: IN_US: %d, IN_byte_0: %u, IN_byte_1: %u, OUT_US: %d, x: %G, y: %G",
                //in_us, in_byte_0, in_byte_1, buf[i], x, y);
    }
    //ALOGE("TESTFILTER frames_num: %u", frames_num);
    //ALOGE("TESTFILTER buffer: %p", buffer);
    //ALOGE("TESTFILTER buf[0]: %p", &(buf[0]));
    //ALOGE("TESTFILTER buf[frames_num]: %p", &(buf[frames_num]));
}

void RecordBufferConverter::reset() {
    if (mResampler != NULL) {
        mResampler->reset();
    }
}

size_t RecordBufferConverter::convert(void *dst,
        AudioBufferProvider *provider, size_t frames)
{
    if (mInputConverterProvider != NULL) {
        mInputConverterProvider->setBufferProvider(provider);
        provider = mInputConverterProvider;
    }

    if (mResampler == NULL) {
        ALOGV("NO RESAMPLING sampleRate:%u mSrcFormat:%#x mDstFormat:%#x",
                mSrcSampleRate, mSrcFormat, mDstFormat);

        AudioBufferProvider::Buffer buffer;
        for (size_t i = frames; i > 0; ) {
            buffer.frameCount = i;
            status_t status = provider->getNextBuffer(&buffer);
            if (status != OK || buffer.frameCount == 0) {
                frames -= i; // cannot fill request.
                break;
            }
            // format convert to destination buffer
            convertNoResampler(dst, buffer.raw, buffer.frameCount);

            dst = (int8_t*)dst + buffer.frameCount * mDstFrameSize;
            i -= buffer.frameCount;
            provider->releaseBuffer(&buffer);
        }
    } else {
         ALOGV("RESAMPLING mSrcSampleRate:%u mDstSampleRate:%u mSrcFormat:%#x mDstFormat:%#x",
                 mSrcSampleRate, mDstSampleRate, mSrcFormat, mDstFormat);

         // reallocate buffer if needed
         if (mBufFrameSize != 0 && mBufFrames < frames) {
             free(mBuf);
             mBufFrames = frames;
             (void)posix_memalign(&mBuf, 32, mBufFrames * mBufFrameSize);
         }
        // resampler accumulates, but we only have one source track
        memset(mBuf, 0, frames * mBufFrameSize);
        frames = mResampler->resample((int32_t*)mBuf, frames, provider);
        // format convert to destination buffer
        convertResampler(dst, mBuf, frames);
    }
    return frames;
}

status_t RecordBufferConverter::updateParameters(
        audio_channel_mask_t srcChannelMask, audio_format_t srcFormat,
        uint32_t srcSampleRate,
        audio_channel_mask_t dstChannelMask, audio_format_t dstFormat,
        uint32_t dstSampleRate)
{
    // quick evaluation if there is any change.
    if (mSrcFormat == srcFormat
            && mSrcChannelMask == srcChannelMask
            && mSrcSampleRate == srcSampleRate
            && mDstFormat == dstFormat
            && mDstChannelMask == dstChannelMask
            && mDstSampleRate == dstSampleRate) {
        return NO_ERROR;
    }

    ALOGV("RecordBufferConverter updateParameters srcMask:%#x dstMask:%#x"
            "  srcFormat:%#x dstFormat:%#x  srcRate:%u dstRate:%u",
            srcChannelMask, dstChannelMask, srcFormat, dstFormat, srcSampleRate, dstSampleRate);
    const bool valid =
            audio_is_input_channel(srcChannelMask)
            && audio_is_input_channel(dstChannelMask)
            && audio_is_valid_format(srcFormat) && audio_is_linear_pcm(srcFormat)
            && audio_is_valid_format(dstFormat) && audio_is_linear_pcm(dstFormat)
            && (srcSampleRate <= dstSampleRate * AUDIO_RESAMPLER_DOWN_RATIO_MAX)
            ; // no upsampling checks for now
    if (!valid) {
        return BAD_VALUE;
    }

    mSrcFormat = srcFormat;
    mSrcChannelMask = srcChannelMask;
    mSrcSampleRate = srcSampleRate;
    mDstFormat = dstFormat;
    mDstChannelMask = dstChannelMask;
    mDstSampleRate = dstSampleRate;

    // compute derived parameters
    mSrcChannelCount = audio_channel_count_from_in_mask(srcChannelMask);
    mDstChannelCount = audio_channel_count_from_in_mask(dstChannelMask);
    mDstFrameSize = mDstChannelCount * audio_bytes_per_sample(mDstFormat);

    // do we need to resample?
    delete mResampler;
    mResampler = NULL;
    if (mSrcSampleRate != mDstSampleRate) {
        mResampler = AudioResampler::create(AUDIO_FORMAT_PCM_FLOAT,
                mSrcChannelCount, mDstSampleRate);
        mResampler->setSampleRate(mSrcSampleRate);
        mResampler->setVolume(AudioMixer::UNITY_GAIN_FLOAT, AudioMixer::UNITY_GAIN_FLOAT);
    }

    // are we running legacy channel conversion modes?
    mIsLegacyDownmix = (mSrcChannelMask == AUDIO_CHANNEL_IN_STEREO
                            || mSrcChannelMask == AUDIO_CHANNEL_IN_FRONT_BACK)
                   && mDstChannelMask == AUDIO_CHANNEL_IN_MONO;
    mIsLegacyUpmix = mSrcChannelMask == AUDIO_CHANNEL_IN_MONO
                   && (mDstChannelMask == AUDIO_CHANNEL_IN_STEREO
                            || mDstChannelMask == AUDIO_CHANNEL_IN_FRONT_BACK);

    // do we need to process in float?
    mRequiresFloat = mResampler != NULL || mIsLegacyDownmix || mIsLegacyUpmix;

    // do we need a staging buffer to convert for destination (we can still optimize this)?
    // we use mBufFrameSize > 0 to indicate both frame size as well as buffer necessity
    if (mResampler != NULL) {
        mBufFrameSize = max(mSrcChannelCount, (uint32_t)FCC_2)
                * audio_bytes_per_sample(AUDIO_FORMAT_PCM_FLOAT);
    } else if (mIsLegacyUpmix || mIsLegacyDownmix) { // legacy modes always float
        mBufFrameSize = mDstChannelCount * audio_bytes_per_sample(AUDIO_FORMAT_PCM_FLOAT);
    } else if (mSrcChannelMask != mDstChannelMask && mDstFormat != mSrcFormat) {
        mBufFrameSize = mDstChannelCount * audio_bytes_per_sample(mSrcFormat);
    } else {
        mBufFrameSize = 0;
    }
    mBufFrames = 0; // force the buffer to be resized.

    // do we need an input converter buffer provider to give us float?
    delete mInputConverterProvider;
    mInputConverterProvider = NULL;
    if (mRequiresFloat && mSrcFormat != AUDIO_FORMAT_PCM_FLOAT) {
        mInputConverterProvider = new ReformatBufferProvider(
                audio_channel_count_from_in_mask(mSrcChannelMask),
                mSrcFormat,
                AUDIO_FORMAT_PCM_FLOAT,
                256 /* provider buffer frame count */);
    }

    // do we need a remixer to do channel mask conversion
    if (!mIsLegacyDownmix && !mIsLegacyUpmix && mSrcChannelMask != mDstChannelMask) {
        (void) memcpy_by_index_array_initialization_from_channel_mask(
                mIdxAry, ARRAY_SIZE(mIdxAry), mDstChannelMask, mSrcChannelMask);
    }
    return NO_ERROR;
}

void RecordBufferConverter::convertNoResampler(
        void *dst, const void *src, size_t frames)
{
    // src is native type unless there is legacy upmix or downmix, whereupon it is float.
    if (mBufFrameSize != 0 && mBufFrames < frames) {
        free(mBuf);
        mBufFrames = frames;
        (void)posix_memalign(&mBuf, 32, mBufFrames * mBufFrameSize);
    }
    // do we need to do legacy upmix and downmix?
    if (mIsLegacyUpmix || mIsLegacyDownmix) {
        void *dstBuf = mBuf != NULL ? mBuf : dst;
        if (mIsLegacyUpmix) {
            upmix_to_stereo_float_from_mono_float((float *)dstBuf,
                    (const float *)src, frames);
        } else /*mIsLegacyDownmix */ {
            downmix_to_mono_float_from_stereo_float((float *)dstBuf,
                    (const float *)src, frames);
        }
        if (mBuf != NULL) {
            memcpy_by_audio_format(dst, mDstFormat, mBuf, AUDIO_FORMAT_PCM_FLOAT,
                    frames * mDstChannelCount);
        }
        return;
    }
    // do we need to do channel mask conversion?
    if (mSrcChannelMask != mDstChannelMask) {
        void *dstBuf = mBuf != NULL ? mBuf : dst;
        memcpy_by_index_array(dstBuf, mDstChannelCount,
                src, mSrcChannelCount, mIdxAry, audio_bytes_per_sample(mSrcFormat), frames);
        if (dstBuf == dst) {
            return; // format is the same
        }
    }
    // convert to destination buffer
    const void *convertBuf = mBuf != NULL ? mBuf : src;
    memcpy_by_audio_format(dst, mDstFormat, convertBuf, mSrcFormat,
            frames * mDstChannelCount);
}

void RecordBufferConverter::convertResampler(
        void *dst, /*not-a-const*/ void *src, size_t frames)
{
    // src buffer format is ALWAYS float when entering this routine
    if (mIsLegacyUpmix) {
        ; // mono to stereo already handled by resampler
    } else if (mIsLegacyDownmix
            || (mSrcChannelMask == mDstChannelMask && mSrcChannelCount == 1)) {
        // the resampler outputs stereo for mono input channel (a feature?)
        // must convert to mono
        downmix_to_mono_float_from_stereo_float((float *)src,
                (const float *)src, frames);
    } else if (mSrcChannelMask != mDstChannelMask) {
        // convert to mono channel again for channel mask conversion (could be skipped
        // with further optimization).
        if (mSrcChannelCount == 1) {
            downmix_to_mono_float_from_stereo_float((float *)src,
                (const float *)src, frames);
        }
        // convert to destination format (in place, OK as float is larger than other types)
        if (mDstFormat != AUDIO_FORMAT_PCM_FLOAT) {
            memcpy_by_audio_format(src, mDstFormat, src, AUDIO_FORMAT_PCM_FLOAT,
                    frames * mSrcChannelCount);
        }
        // channel convert and save to dst
        memcpy_by_index_array(dst, mDstChannelCount,
                src, mSrcChannelCount, mIdxAry, audio_bytes_per_sample(mDstFormat), frames);
        return;
    }
    // convert to destination format and save to dst
    memcpy_by_audio_format(dst, mDstFormat, src, AUDIO_FORMAT_PCM_FLOAT,
            frames * mDstChannelCount);
}

// ----------------------------------------------------------------------------
} // namespace android
