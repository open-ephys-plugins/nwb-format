/*
 ------------------------------------------------------------------

 This file is part of the Open Ephys GUI
 Copyright (C) 2014 Open Ephys

 ------------------------------------------------------------------

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */
 
 #include "NWBRecording.h"
#define MAX_BUFFER_SIZE 40960
 
 using namespace NWBRecording;
 
 NWBRecordEngine::NWBRecordEngine() 
 {
	 smpBuffer.malloc(MAX_BUFFER_SIZE);
     
     
 }
 
 NWBRecordEngine::~NWBRecordEngine()
 {
 }
 
 String NWBRecordEngine::getEngineId() const
 {
	 return "NWB2"; //a text identifier
 }
 
 void NWBRecordEngine::openFiles(File rootFolder, int experimentNumber, int recordingNumber)
 {
	 
	 //Called when acquisition starts, to open the files
	 String basepath = rootFolder.getFullPathName() + rootFolder.getSeparatorString() + "experiment" + String(experimentNumber) + ".nwb";
     
     if (!File(basepath).exists())
     {
         Uuid identifier;
         identifierText = identifier.toString();
     }
     
	 recordFile = new NWBFile(basepath, CoreServices::getGUIVersion(), identifierText);
	 recordFile->setXmlText(getLatestSettingsXml());

	 datasetIndexes.insertMultiple(0, 0, getNumRecordedContinuousChannels());
	 writeChannelIndexes.insertMultiple(0, 0, getNumRecordedContinuousChannels());
     continuousChannelGroups.clear();

	 //Generate the continuous datasets info array, seeking for different combinations of recorded processor and source processor
     Array<const ContinuousChannel*> firstChannels;
     Array<int> channelCounts;
     
     int streamIndex = -1;
     int streamChannelCount = 0;
	 
     for (int ch = 0; ch < getNumRecordedContinuousChannels(); ch++)
     {

         int globalIndex = getGlobalIndex(ch); // the global channel index
         int localIndex = getLocalIndex(ch); // the local channel index (within a stream)
         
         const ContinuousChannel* channelInfo = getContinuousChannel(globalIndex); // channel info object
		 
         int sourceId = channelInfo->getSourceNodeId();
         int streamId = channelInfo->getStreamId();

         if (localIndex == 0)
         {
             firstChannels.add(channelInfo);
             streamIndex++;
             
             if (streamIndex > 0)
             {
                 channelCounts.add(streamChannelCount);
             }
                 
             ContinuousGroup newGroup;
             newGroup.add(channelInfo);
             continuousChannelGroups.add(newGroup);

             streamChannelCount = 0;
         }
         
         datasetIndexes.set(ch, streamIndex);
         writeChannelIndexes.set(ch, localIndex);
	 }
     
	 int nEvents = getNumRecordedEventChannels();
	 for (int i = 0; i < nEvents; i++)
		 eventChannels.add(getEventChannel(i));

	 int nSpikes = getNumRecordedSpikeChannels();
	 for (int i = 0; i < nSpikes; i++)
		 spikeChannels.add(getSpikeChannel(i));

	 //open the file
	 recordFile->open(getNumRecordedContinuousChannels() + continuousChannelGroups.size() + eventChannels.size() + spikeChannels.size()); //total channels + timestamp arrays, to create a big enough buffer

	 //create the recording
	 recordFile->startNewRecording(recordingNumber, continuousChannelGroups, eventChannels, spikeChannels);
	
 }

 
 void NWBRecordEngine::closeFiles()
 {
	 //Called when acquisition stops. Should close the files and leave the processor in a reset status
	 recordFile->stopRecording();
	 recordFile->close();
	 recordFile = nullptr;
     spikeChannels.clear();
     eventChannels.clear();
     continuousChannelGroups.clear();
     datasetIndexes.clear();
     writeChannelIndexes.clear();
 }

 

void NWBRecordEngine::writeContinuousData(int writeChannel,
                                          int realChannel,
                                          const float* dataBuffer,
                                          const double* timestampBuffer,
                                          int size)
{
    recordFile->writeData(datasetIndexes[writeChannel], writeChannelIndexes[writeChannel], size, dataBuffer, getContinuousChannel(realChannel)->getBitVolts());

    /* All channels in a dataset have the same number of samples and share timestamps. But since this method is called asynchronously, the timestamps might not be in sync during acquisition, so we chose a channel and write the timestamps when writing that channel's data */
    if (writeChannelIndexes[writeChannel] == 0)
    {
        int64 baseTS = getTimestamp(writeChannel);
        
        for (int i = 0; i < size; i++)
        {
            smpBuffer[i] = baseTS + i;
        }
        
        recordFile->writeTimestamps(datasetIndexes[writeChannel], size, timestampBuffer);
        recordFile->writeSampleNumbers(datasetIndexes[writeChannel], size, smpBuffer);
    }
}
 
void NWBRecordEngine::writeEvent(int eventIndex, const MidiMessage& event) 
{
	const EventChannel* channel = getEventChannel(eventIndex);
	EventPtr eventStruct = Event::deserialize(event, channel);

	recordFile->writeEvent(eventIndex, channel, eventStruct);
}

void NWBRecordEngine::writeTimestampSyncText(uint64 streamId, int64 timestamp, float sourceSampleRate, String text)
{
	recordFile->writeTimestampSyncText(streamId, timestamp, sourceSampleRate, text);
}


void NWBRecordEngine::writeSpike(int electrodeIndex, const Spike* spike) 
{
	const SpikeChannel* channel = getSpikeChannel(electrodeIndex);

	recordFile->writeSpike(electrodeIndex, channel, spike);
}

RecordEngineManager* NWBRecordEngine::getEngineManager()
{
	//static factory that instantiates the engine manager, which allows to configure recording options among other things. See OriginalRecording to see how to create options for a record engine
	RecordEngineManager* man = new RecordEngineManager("NWB2", "NWB2", &(engineFactory<NWBRecordEngine>));
	EngineParameter* param;
	param = new EngineParameter(EngineParameter::STR, 0, "Identifier Text", String());
	man->addParameter(param);
	return man;
	
}

void NWBRecordEngine::setParameter(EngineParameter& parameter)
{
	strParameter(0, identifierText);
}
