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

#ifndef NWBFORMAT_H
#define NWBFORMAT_H

#include <HDF5FileFormat.h>
#include <RecordingLib.h>
#include <ProcessorHeaders.h>

using namespace OpenEphysHDF5;

namespace NWBRecording
{

    typedef Array<const ContinuousChannel*> ContinuousGroup;

    /**
        Represents an NWB TimeSeries dataset
     */
	class TimeSeries
	{
	public:
        
        /** Holds the sample data */
		ScopedPointer<HDF5RecordingData> baseDataSet;
        
        /** Holds the timestamps (in seconds) for each sample */
		ScopedPointer<HDF5RecordingData> timestampDataSet;
        
        /** Holds the sample number for each sample (relative to the start of acquisition) */
		ScopedPointer<HDF5RecordingData> sampleNumberDataSet;
        
        /** Holds the electrode index for each channel */
		ScopedPointer<HDF5RecordingData> electrodeDataSet;
        
        /** Holds the TTL word for each sample (TTL event TimeSeries only) */
		ScopedPointer<HDF5RecordingData> ttlWordDataSet;
        
        /** Holds the metadata for each event (if applicable) */
		OwnedArray<HDF5RecordingData> metaDataSet;
        
        /** The path to this dataset within the NWB file */
		String basePath;
        
        /** Total number of samples written */
		uint64 numSamples{ 0 };
	};

    /**
        
        Represents an NWB 2.0 File (a specific type of HDF5 file)
            
     */
	class NWBFile : public HDF5FileBase
	{
	public:
        
        /** Constructor */
		NWBFile(String fName, String ver, String identifier);
        
        /** Destructor */
        ~NWBFile() { }
        
        /** Creates the groups required for a new recording, given an array of continuous channels, event channels, and spike channels*/
		bool startNewRecording(int recordingNumber,
                               const Array<ContinuousGroup>& continuousArray,
                               const Array<const EventChannel*>& eventArray,
                               const Array<const SpikeChannel*>& electrodeArray);
        
        /** Writes the num_samples value and closes the relevent datasets */
		void stopRecording();
        
        /** Writes continuous data for a particular channel */
		void writeData(int datasetID, int channel, int nSamples, const float* data, float bitVolts);
        
        /** Writes synchronized timestamps for a particular continuous dataset */
		void writeTimestamps(int datasetID, int nSamples, const double* data);
        
        /** Writes sample numbers for a particular continuous dataset */
		void writeSampleNumbers(int datasetID, int nSamples, const int64* data);
        
        /** Writes electrode numbers for a continuous dataset */
		void writeElectrodes(int datasetID, int start, int nElectrodes);
        
        /** Writes a spike event*/
		void writeSpike(int electrodeId, const SpikeChannel* channel, const Spike* event);
        
        /** Writes an event (TEXT or TTL) */
		void writeEvent(int eventID, const EventChannel* channel, const Event* event);
        
        /** Writes a timestamp sync text event */
		void writeTimestampSyncText(uint16 sourceID,
                                    int64 timestamp,
                                    float sourceSampleRate,
                                    String text);
        
        /** Returns the name of this NWB file */
		String getFileName() override;

	protected:
        
        /** Initializes the default groups */
		int createFileStructure() override;

	private:

        /** Creates a new dataset to hold text data (messages) */
		void createTextDataSet(String path, String name, String text);
        
        /** Creates a new dataset to hold binary events */
		void createBinaryDataSet(String path, String name, HDF5FileBase::BaseDataType type, int length, void* data);
        
        /** Returns the HDF5 data type for a given event channel type */
		static HDF5FileBase::BaseDataType getEventH5Type(EventChannel::Type type, int length = 1);
		
        /** Returns the HDF5 data type for a given metadata type*/
        static HDF5FileBase::BaseDataType getMetadataH5Type(MetadataDescriptor::MetadataType type, int length = 1);

        /** Creates a time series dataset*/
		bool createTimeSeriesBase(String basePath, String description, String neurodata_type);
		
        /** Creates dataset attributes */
        bool createExtraInfo(String basePath, String name, String desc, String id, uint16 index, uint16 typeIndex);
        
        /** Creates a dataset of synchronized timestamps */
		HDF5RecordingData* createTimestampDataSet(String basePath, int chunk_size);
        
        /** Creates a dataset of sample numbers */
		HDF5RecordingData* createSampleNumberDataSet(String basePath, int chunk_size);
        
        /** Creates a dataset for electrode indices */
		HDF5RecordingData* createElectrodeDataSet(String basePath, String description, int chunk_size);
		
        /** Adds attributes (e.g. conversion, resolution) to a continuous dataset */
        void createDataAttributes(String basePath, float conversion, float resolution, String unit);
		
        /** Creates a dataset for channel metdata */
        bool createChannelMetadataSets(String basePath, const MetadataObject* info);
		
        /** Creates a dataset for event metdata */
        bool createEventMetadataSets(String basePath, TimeSeries* timeSeries, const MetadataEventObject* info);

        /** Writes metadata associated with an event*/
		void writeEventMetadata(TimeSeries* timeSeries, const MetadataEventObject* info, const MetadataEvent* event);
		
		const String filename;
		const String GUIVersion;

		OwnedArray<TimeSeries>  continuousDataSets;
		OwnedArray<TimeSeries> spikeDataSets;
		OwnedArray<TimeSeries> eventDataSets;
		ScopedPointer<TimeSeries> syncMsgDataSet;

		const String identifierText;

		HeapBlock<float> scaledBuffer;
		HeapBlock<int16> intBuffer;
		size_t bufferSize;

		JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(NWBFile);

	};

}

#endif
