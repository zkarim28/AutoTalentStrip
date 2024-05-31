#include "PluginProcessor.h"
#include "PluginEditor.h"

AutoTalentStripPluginAudioProcessorEditor::AutoTalentStripPluginAudioProcessorEditor (AutoTalentStripPluginAudioProcessor& p)
    : AudioProcessorEditor (&p), processor (p)
{
    setSize (400, 300);
}

AutoTalentStripPluginAudioProcessorEditor::~AutoTalentStripPluginAudioProcessorEditor()
{
}

void AutoTalentStripPluginAudioProcessorEditor::paint (juce::Graphics& g)
{
    g.fillAll (getLookAndFeel().findColour (juce::ResizableWindow::backgroundColourId));
}

void AutoTalentStripPluginAudioProcessorEditor::resized()
{
}
