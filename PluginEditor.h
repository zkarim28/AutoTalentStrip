#pragma once

#include <JuceHeader.h>
#include "PluginProcessor.h"

class AutoTalentStripPluginAudioProcessorEditor  : public juce::AudioProcessorEditor
{
public:
    AutoTalentStripPluginAudioProcessorEditor (AutoTalentStripPluginAudioProcessor&);
    ~AutoTalentStripPluginAudioProcessorEditor() override;

    void paint (juce::Graphics&) override;
    void resized() override;

private:
    AutoTalentStripPluginAudioProcessor& processor;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (AutoTalentStripPluginAudioProcessorEditor)
};
