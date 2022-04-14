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
 
 #include "NWBFormat.h"
 
 using namespace NWBRecording;

#ifndef EVENT_CHUNK_SIZE
#define EVENT_CHUNK_SIZE 8
#endif

#ifndef SPIKE_CHUNK_XSIZE
#define SPIKE_CHUNK_XSIZE 8
#endif

#ifndef SPIKE_CHUNK_YSIZE
#define SPIKE_CHUNK_YSIZE 40
#endif

 #define MAX_BUFFER_SIZE 40960

NWBFile::NWBFile(String fName, String ver, String idText) :
    HDF5FileBase(),
    filename(fName),
    identifierText(idText),
    GUIVersion(ver)
{
	 readyToOpen = true; //In KWIK this is in initFile, but the new recordEngine methods make it safe for it to be here

	 scaledBuffer.malloc(MAX_BUFFER_SIZE);
	 intBuffer.malloc(MAX_BUFFER_SIZE);
	 bufferSize = MAX_BUFFER_SIZE;
}

NWBFile::~NWBFile()
{
	continuousDataSets.clear();
	spikeDataSets.clear();
	eventDataSets.clear();
	syncMsgDataSet.reset();
}
 
int NWBFile::createFileStructure()
{

	setAttributeStr("core", "/", "namespace");
	setAttributeStr("NWBFile", "/", "neurodata_type");
	setAttributeStr("2.4.0", "/", "nwb_version");
	setAttributeStr(identifierText, "/", "object_id");

	if (createGroup("/acquisition")) return -1;

	if (createGroup("/analysis")) return -1;

	String time = Time::getCurrentTime().formatted("%Y-%m-%dT%H:%M:%S") + Time::getCurrentTime().getUTCOffsetString(true);

	createTextDataSet("", "file_create_date", time);

	if (createGroup("/general")) return -1;
	if (createGroup("general/devices")) return -1;
	if (createGroup("general/extracellular_ephys")) return -1;
	if (createGroup("/processing")) return -1;

	if (createGroup("/stimulus")) return -1;
	if (createGroup("/stimulus/presentation")) return -1;
	if (createGroup("/stimulus/templates")) return -1;

	createStringDataSet("/session_description", "Recording with the Open Ephys GUI");
	createStringDataSet("/session_start_time", time);
	createStringDataSet("/timestamps_reference_time", time);
	createStringDataSet("/identifier", "test-identifier");

	return 0;

}
 
bool NWBFile::startNewRecording(
	int recordingNumber, 
	const Array<ContinuousGroup>& continuousArray,
	const Array<const EventChannel*>& eventArray, 
	const Array<const SpikeChannel*>& electrodeArray)
{

    // all recorded data is stored in the "acquisition" group
	String rootPath = "/acquisition/";
	
    // sub-path within "acquisition" group
	String basePath;

	continuousDataSets.clearQuick(true);
	spikeDataSets.clearQuick(true);
	eventDataSets.clearQuick(true);

	ScopedPointer<TimeSeries> timeSeries;
	ScopedPointer<HDF5RecordingData> dSet;
	ScopedPointer<HDF5RecordingData> eSet;
	ScopedPointer<HDF5RecordingData> starting_time;

    int nRecordedStreams = continuousArray.size();
	int totalElectrodeCount = 0;

	for (int i = 0; i < nRecordedStreams; i++)
	{

		// All channels in a group will share the same source information 
		//  (any caller to this method MUST assure this happen), so we just need to look at the first channel
		const ContinuousChannel* info = continuousArray.getReference(i)[0];

		String desc = info->getSourceNodeName() + "-" 
					+ String(info->getSourceNodeId())
					+  "." + info->getStreamName();
        
		basePath = rootPath + desc;
        
        if (recordingNumber == 0)
            if (!createTimeSeriesBase(basePath, "Stores voltage data from extracellular recordings", "")) return false;

        timeSeries = new TimeSeries();
        timeSeries->basePath = basePath;
		
		//std::cout << basePath << "/data" << std::endl;
        
        dSet = createDataSet(BaseDataType::I16,
                                0,
                                continuousArray.getReference(i).size(),
                                CHUNK_XSIZE,
                                basePath + "/data");
            
        if (dSet == nullptr)
        {
            std::cerr << "Error creating dataset for " << desc << std::endl;
            return false;
        }
        else
        {
            createDataAttributes(basePath, info->getBitVolts(), info->getBitVolts() / 65536, info->getUnits());
        }

		
        timeSeries->baseDataSet = dSet;

        dSet = createTimestampDataSet(basePath + "/timestamps", CHUNK_XSIZE);
        if (dSet == nullptr) return false;
        
        timeSeries->timestampDataSet = dSet;
        
		//std::cout << basePath << "/sample_numbers" << std::endl;
       
        dSet = createSampleNumberDataSet(basePath + "/sample_numbers", CHUNK_XSIZE);
        if (dSet == nullptr) return false;

        timeSeries->sampleNumberDataSet = dSet;
        
		//std::cout << basePath << "/electrodes" << std::endl;
       
        int numElectrodesInStream = continuousArray.getReference(i).size();
        
        dSet = createElectrodeDataSet(basePath + "/electrodes", desc, CHUNK_XSIZE);
        if (dSet == nullptr) return false;
        writeElectrodes(i, totalElectrodeCount, numElectrodesInStream);
        
        totalElectrodeCount += numElectrodesInStream;
		
		timeSeries->electrodeDataSet = dSet;

		continuousDataSets.add(timeSeries.release());
	}

	//std::cout << "Created continuous channels " << std::endl;

	int nRecordedSpikeElectrodes;
	nRecordedSpikeElectrodes = electrodeArray.size();
	std::unordered_set<String> spikeStreams;

	String currentGroup = "";
	
	for (int i = 0; i < nRecordedSpikeElectrodes; i++)
	{
		const SpikeChannel* sourceInfo = electrodeArray[i];

		String sourceName = sourceInfo->getSourceNodeName() + "-" + String(sourceInfo->getSourceNodeId());
		sourceName += "." + sourceInfo->getStreamName();
        
        //if (recordingNumber > 0)
        //    sourceName += "." + String(recordingNumber + 1);
        
		basePath = rootPath + sourceName + ".spikes";

        createGroupIfDoesNotExist(basePath);

		basePath += "/" + sourceInfo->getName();
        
        createGroupIfDoesNotExist(basePath);

        if (!createTimeSeriesBase(basePath, "Stores acquired spike data from extracellular recordings", "")) return false;

		timeSeries = new TimeSeries();
        timeSeries->basePath = basePath;
        
        dSet = createDataSet(BaseDataType::I16, 0, sourceInfo->getNumChannels(), sourceInfo->getTotalSamples(), SPIKE_CHUNK_XSIZE, basePath + "/data");
        if (dSet == nullptr)
        {
            std::cerr << "Error creating dataset for electrode " << i << std::endl;
            return false;
        }
        else
        {
            createDataAttributes(basePath, sourceInfo->getChannelBitVolts(0), sourceInfo->getChannelBitVolts(0) / 65536, "volt");
        }

        timeSeries->baseDataSet = dSet;
       
        dSet = createTimestampDataSet(basePath + "/timestamps", SPIKE_CHUNK_XSIZE);
		if (dSet == nullptr) return false;
        
        timeSeries->timestampDataSet = dSet;
    
		spikeDataSets.add(timeSeries.release());

	}

	//std::cout << "Created spike channels " << std::endl;
	
	int nEvents = eventArray.size();
	int nTTL = 0;
	int nTXT = 0;
	int nBIN = 0;
	
	for (int i = 0; i < nEvents; i++)
	{

		basePath = rootPath;
		const EventChannel* info = eventArray[i];
		String sourceName = info->getSourceNodeName() + "-" + String(info->getSourceNodeId());
		sourceName += "." + info->getStreamName();
        
        //if (recordingNumber > 0)
         //   sourceName += "." + String(recordingNumber + 1);

		String series;

		String helpText;

		switch (info->getType())
		{	
		case EventChannel::TTL:
			nTTL += 1;
			basePath += sourceName + ".TTL";
			series = "IntervalSeries";
			helpText = "Stores the start and stop times for TTL events";
			break;
		case EventChannel::TEXT:
			nTXT += 1;
			basePath += "messages";
			//if (recordingNumber > 0)
			//	basePath += "." + String(recordingNumber + 1);
			series = "AnnotationSeries";
			helpText = "Time-stamped annotations about an experiment";
			break;
		default:
			nBIN += 1;
			basePath += sourceName + ".custom";
			series = "IntervalSeries";
			helpText = "Stores arbitrary binary data";
			break;
		 }

        if (recordingNumber == 0)
            if (!createTimeSeriesBase(basePath, helpText, series)) return false;

		timeSeries = new TimeSeries();
        timeSeries->basePath = basePath;
        
        if (info->getType() >= EventChannel::BinaryDataType::BINARY_BASE_VALUE) //only binary events have length greater than 1
        {
            dSet = createDataSet(getEventH5Type(info->getType(), info->getLength()), 0, info->getLength(), EVENT_CHUNK_SIZE, basePath + "/data");;
        }
        else
        {
            dSet = createDataSet(getEventH5Type(info->getType(), info->getLength()), 0, EVENT_CHUNK_SIZE, basePath + "/data");
        }

        if (dSet == nullptr)
        {
            std::cerr << "Error creating dataset for event " << info->getName() << std::endl;
            return false;
        }
        else
        {
            createDataAttributes(basePath, NAN, NAN, "n/a");
        }


        timeSeries->baseDataSet = dSet;

        dSet = createTimestampDataSet(basePath + "/timestamps", EVENT_CHUNK_SIZE);
        if (dSet == nullptr) return false;

        timeSeries->timestampDataSet = dSet;

		dSet = createSampleNumberDataSet(basePath + "/sample_numbers", EVENT_CHUNK_SIZE);
		if (dSet == nullptr) return false;

		timeSeries->sampleNumberDataSet = dSet;

		if (info->getType() == EventChannel::TTL)
		{
            
            dSet = createDataSet(BaseDataType::U64, 0, info->getDataSize(), EVENT_CHUNK_SIZE, basePath + "/full_word");
            if (dSet == nullptr) return false;
			
            timeSeries->ttlWordDataSet = dSet;
		}
        
		eventDataSets.add(timeSeries.release());

	}

	//std::cout << "Created event channels " << std::endl;
	
	basePath = rootPath + "sync_messages";

	String desc = "Stores recording start timestamps for each processor in text format";
	
    syncMsgDataSet = std::make_unique<TimeSeries>();
	syncMsgDataSet->basePath = basePath;

    if (recordingNumber == 0)
    {

		if (!createTimeSeriesBase(basePath, "Auto-generated messages at the start of each recording", "AnnotationSeries")) return false;

        dSet = createDataSet(BaseDataType::STR(100), 0, 1, basePath + "/data");
  
        if (dSet == nullptr)
        {
            std::cerr << "Error creating dataset for sync messages" << std::endl;
            return false;
        }
        else
        {
            createDataAttributes(basePath, NAN, NAN, "n/a");
        }
	}
	else {
		dSet = getDataSet(basePath + "/data");
	}
   
	syncMsgDataSet->baseDataSet = dSet;

    if (recordingNumber == 0)
    {
        dSet = createSampleNumberDataSet(basePath + "/sample_numbers", 1);
        if (dSet == nullptr) return false;
	}
	else {
		dSet = getDataSet(basePath + "/sample_numbers");
	}
    
	syncMsgDataSet->sampleNumberDataSet = dSet;

	if (recordingNumber == 0)
	{
		dSet = createTimestampDataSet(basePath + "/timestamps", 1);
		if (dSet == nullptr) return false;
	}
	else {
		dSet = getDataSet(basePath + "/timestamps");
	}

	syncMsgDataSet->timestampDataSet = dSet;

	return true;

 }
 
 void NWBFile::stopRecording()
 {

	 const TimeSeries* tsStruct;

	 for (int i = 0; i < continuousDataSets.size(); i++)
	 {
		 tsStruct = continuousDataSets[i];
		 CHECK_ERROR(setAttribute(BaseDataType::U64, &(tsStruct->numSamples), tsStruct->basePath, "num_samples"));
	 }

	 for (int i = 0; i < spikeDataSets.size(); i++)
	 {
		 tsStruct = spikeDataSets[i];
		 CHECK_ERROR(setAttribute(BaseDataType::U64, &(tsStruct->numSamples), tsStruct->basePath, "num_samples"));
	 }
	 
	 for (int i = 0; i < eventDataSets.size(); i++)
	 {
		 tsStruct = eventDataSets[i];
		 CHECK_ERROR(setAttribute(BaseDataType::U64, &(tsStruct->numSamples), tsStruct->basePath, "num_samples"));
	 }
	 
	 CHECK_ERROR(setAttribute(BaseDataType::U64, &(syncMsgDataSet->numSamples), syncMsgDataSet->basePath, "num_samples"));
 }
 
 void NWBFile::writeData(int datasetID, int channel, int nSamples, const float* data, float bitVolts)
 {
	 if (!continuousDataSets[datasetID])
		 return;

	 if (nSamples > bufferSize) //Shouldn't happen, and if it happens it'll be slow, but better this than crashing. Will be reset on file close and reset.
	 {
		 std::cerr << "Write buffer overrun, resizing to" << nSamples << std::endl;
		 bufferSize = nSamples;
		 scaledBuffer.malloc(nSamples);
		 intBuffer.malloc(nSamples);
	 }

	 double multFactor = 1 / (float(0x7fff) * bitVolts);
	 FloatVectorOperations::copyWithMultiply(scaledBuffer.getData(), data, multFactor, nSamples);
	 AudioDataConverters::convertFloatToInt16LE(scaledBuffer.getData(), intBuffer.getData(), nSamples);

	 continuousDataSets[datasetID]->baseDataSet->writeDataRow(channel, nSamples, BaseDataType::I16, intBuffer);
	 //CHECK_ERROR();
	 
	 /* Since channels are filled asynchronouysly by the Record Thread, there is no guarantee
		that at a any point in time all channels in a dataset have the same number of filled samples.
		However, since each dataset is filled from a single source, all channels must have the
		same number of samples at acquisition stop. To keep track of the written samples we must chose
		an arbitrary channel, and at the end all channels will be the same. */

	 if (channel == 0) //there will always be a first channel or there wouldn't be dataset
		 continuousDataSets[datasetID]->numSamples += nSamples;
 }

 void NWBFile::writeSampleNumbers(int datasetID, int nSamples, const int64* data)
 {
	 if (!continuousDataSets[datasetID])
		 return;

	 CHECK_ERROR(continuousDataSets[datasetID]->sampleNumberDataSet->writeDataBlock(nSamples, BaseDataType::I64, data));
 }

 void NWBFile::writeTimestamps(int datasetID, int nSamples, const double* data)
 {
	 if (!continuousDataSets[datasetID])
		 return;

	 CHECK_ERROR(continuousDataSets[datasetID]->timestampDataSet->writeDataBlock(nSamples, BaseDataType::F64, data));
 }

 void NWBFile::writeElectrodes(int datasetID, int start, int nElectrodes)
 {
	 if (!continuousDataSets[datasetID])
		 return;

	std::vector<int> electrodeNumbers;
	for (int i = start; i < start + nElectrodes; i++)
		electrodeNumbers.push_back(i);

	 CHECK_ERROR(continuousDataSets[datasetID]->electrodeDataSet->writeDataBlock(nElectrodes, BaseDataType::I32, &electrodeNumbers[0]));
 }

 void NWBFile::writeSpike(int electrodeId, const SpikeChannel* channel, const Spike* event)
 {
	 if (!spikeDataSets[electrodeId])
		 return;
	 int nSamples = channel->getTotalSamples() * channel->getNumChannels();

	 if (nSamples > bufferSize) //Shouldn't happen, and if it happens it'll be slow, but better this than crashing. Will be reset on file close and reset.
	 {
		 std::cerr << "Write buffer overrun, resizing to" << nSamples << std::endl;
		 bufferSize = nSamples;
		 scaledBuffer.malloc(nSamples);
		 intBuffer.malloc(nSamples);
	 }

	 double multFactor = 1 / (float(0x7fff) * channel->getChannelBitVolts(0));
	 FloatVectorOperations::copyWithMultiply(scaledBuffer.getData(), event->getDataPointer(), multFactor, nSamples);
	 AudioDataConverters::convertFloatToInt16LE(scaledBuffer.getData(), intBuffer.getData(), nSamples);

	 double timestampSec = event->getTimestampInSeconds();

	 CHECK_ERROR(spikeDataSets[electrodeId]->baseDataSet->writeDataBlock(1, BaseDataType::I16, intBuffer));
	 CHECK_ERROR(spikeDataSets[electrodeId]->timestampDataSet->writeDataBlock(1, BaseDataType::F64, &timestampSec));
	 writeEventMetadata(spikeDataSets[electrodeId], channel, event);

	 spikeDataSets[electrodeId]->numSamples += 1;

 }

 void NWBFile::writeEvent(int eventID, const EventChannel* channel, const Event* event)
 {
	 if (!eventDataSets[eventID])
		 return;
	 
	 const void* dataSrc;
	 BaseDataType type;
	 int8 ttlVal;
	 String text;

	 switch (event->getEventType())
	 {
	 case EventChannel::TTL:
		 ttlVal = (static_cast<const TTLEvent*>(event)->getState() ? 1 : -1) * (static_cast<const TTLEvent*>(event)->getLine() + 1);
		 dataSrc = &ttlVal;
		 type = BaseDataType::I8;
		 break;
	 case EventChannel::TEXT:
		 text = static_cast<const TextEvent*>(event)->getText();
		 dataSrc = text.toUTF8().getAddress();
		 type = BaseDataType::STR(text.length());
		 break;
	 default:
		 dataSrc = static_cast<const BinaryEvent*>(event)->getBinaryDataPointer();
		 type = getEventH5Type(event->getEventType());
		 break;
	 }
	 CHECK_ERROR(eventDataSets[eventID]->baseDataSet->writeDataBlock(1, type, dataSrc));

	 const double timeSec = event->getTimestampInSeconds();

	 CHECK_ERROR(eventDataSets[eventID]->timestampDataSet->writeDataBlock(1, BaseDataType::F64, &timeSec));

	 const int64 sampleNumber = event->getSampleNumber();

	 CHECK_ERROR(eventDataSets[eventID]->sampleNumberDataSet->writeDataBlock(1, BaseDataType::I64, &sampleNumber));

	 if (event->getEventType() == EventChannel::TTL)
	 {
        const uint64 ttlWord = static_cast<const TTLEvent*>(event)->getWord();
		CHECK_ERROR(eventDataSets[eventID]->ttlWordDataSet->writeDataBlock(1, BaseDataType::U64, &ttlWord));
	 }
	 
	 eventDataSets[eventID]->numSamples += 1;
 }

 void NWBFile::writeTimestampSyncText(uint16 sourceID, int64 sampleNumber, float sourceSampleRate, String text)
 {
	 CHECK_ERROR(syncMsgDataSet->baseDataSet->writeDataBlock(1, BaseDataType::STR(text.length()), text.toUTF8()));
     
	 CHECK_ERROR(syncMsgDataSet->sampleNumberDataSet->writeDataBlock(1, BaseDataType::I64, &sampleNumber));

	 double timestamp = (double)sampleNumber;

	 CHECK_ERROR(syncMsgDataSet->timestampDataSet->writeDataBlock(1, BaseDataType::F64, &timestamp));

	 syncMsgDataSet->numSamples += 1;
 }

 
 String NWBFile::getFileName()
 {
	 return filename;
 }

  bool NWBFile::createTimeSeriesBase(String basePath, String description, String neurodata_type)
 {
	 if (createGroup(basePath)) return false;
	 CHECK_ERROR(setAttributeStr(" ", basePath, "comments"));
	 CHECK_ERROR(setAttributeStr(description, basePath, "description"));
	 CHECK_ERROR(setAttributeStr("core", basePath, "namespace"));
	 CHECK_ERROR(setAttributeStr(neurodata_type, basePath, "neurodata_type"));
	 return true;
 }

  void NWBFile::createDataAttributes(String basePath, float conversion, float resolution, String unit)
  {
		  CHECK_ERROR(setAttribute(BaseDataType::F32, &conversion, basePath + "/data", "conversion"));
		  CHECK_ERROR(setAttribute(BaseDataType::F32, &resolution, basePath + "/data", "resolution"));
		  CHECK_ERROR(setAttributeStr(unit, basePath + "/data", "unit"));
  }

  HDF5RecordingData* NWBFile::createTimestampDataSet(String path, int chunk_size)
  {
	  HDF5RecordingData* tsSet = createDataSet(BaseDataType::F64, 0, chunk_size, path);
	  if (!tsSet)
		  std::cerr << "Error creating timestamp dataset in " << path << std::endl;
	  else
	  {
		  const int32 one = 1;
		  CHECK_ERROR(setAttribute(BaseDataType::I32, &one, path, "interval"));
		  CHECK_ERROR(setAttributeStr("seconds", path, "unit"));
	  }
	  return tsSet;
  }

   HDF5RecordingData* NWBFile::createSampleNumberDataSet(String path, int chunk_size)
  {
	  HDF5RecordingData* tsSet = createDataSet(BaseDataType::I64, 0, chunk_size, path);
	  if (!tsSet)
		  std::cerr << "Error creating sample number dataset in " << path << std::endl;
	  else
	  {
		  const int32 one = 1;
		  CHECK_ERROR(setAttribute(BaseDataType::I32, &one, path, "interval"));
		  CHECK_ERROR(setAttributeStr("samples", path, "unit"));
	  }
	  return tsSet;
  }

  HDF5RecordingData *NWBFile::createElectrodeDataSet(String path, String description, int chunk_size)
{
	HDF5RecordingData *elSet = createDataSet(BaseDataType::I32, 1, chunk_size, path);
	if (!elSet)
		std::cerr << "Error creating electrode dataset in " << path << std::endl;
	else
	{
		CHECK_ERROR(setAttributeStr(description, path, "description"));
		CHECK_ERROR(setAttributeStr("hdmf-common", path, "namespace"));
		CHECK_ERROR(setAttributeStr("DynamicTableRegion", path, "neurodata_type"));
		//TODO: object_id(uuid), table
	}
	return elSet;
}

  bool NWBFile::createExtraInfo(String basePath, String name, String desc, String id, uint16 index, uint16 typeIndex)
  {
	  if (createGroup(basePath)) return false;
	  CHECK_ERROR(setAttributeStr("openephys:<channel_info>/", basePath, "schema_id"));
	  CHECK_ERROR(setAttributeStr(name, basePath, "name"));
	  CHECK_ERROR(setAttributeStr(desc, basePath, "description"));
	  CHECK_ERROR(setAttributeStr(id, basePath, "identifier"));
	  CHECK_ERROR(setAttribute(BaseDataType::U16, &index, basePath, "source_index"));
	  CHECK_ERROR(setAttribute(BaseDataType::U16, &typeIndex, basePath, "source_type_index"));
	  return true;
  }

  bool NWBFile::createChannelMetadataSets(String basePath, const MetadataObject* info)
  {
	  if (!info) return false;
	  if (createGroup(basePath)) return false;
	  CHECK_ERROR(setAttributeStr("openephys:<metadata>/", basePath, "schema_id"));
	  int nMetadata = info->getMetadataCount();
	  
	  for (int i 
	  	= 0; i < nMetadata; i++)
	  {
		  const MetadataDescriptor* desc = info->getMetadataDescriptor(i);
		  String fieldName = "Field_" + String(i+1);
		  String name = desc->getName();
		  String description = desc->getDescription();
		  String identifier = desc->getIdentifier();
		  BaseDataType type = getMetadataH5Type(desc->getType(), desc->getLength()); //only string types use length, for others is always set to 1. If array types are implemented, change this
		  int length = desc->getType() == MetadataDescriptor::CHAR ? 1 : desc->getLength(); //strings are a single element of length set in the type (see above) while other elements are saved a
		  HeapBlock<char> data(desc->getDataSize());
		  info->getMetadataValue(i)->getValue(static_cast<void*>(data.getData()));
		  createBinaryDataSet(basePath, fieldName, type, length, data.getData());
		  String fullPath = basePath + "/" + fieldName;
		  CHECK_ERROR(setAttributeStr("openephys:<metadata>/", fullPath, "schema_id"));
		  CHECK_ERROR(setAttributeStr(name, fullPath, "name"));
		  CHECK_ERROR(setAttributeStr(description, fullPath, "description"));
		  CHECK_ERROR(setAttributeStr(identifier, fullPath, "identifier"));
	  }
	  return true;
  }

 
  bool NWBFile::createEventMetadataSets(String basePath, TimeSeries* timeSeries, const MetadataEventObject* info)
  {
	  if (!info) return false;
	  if (createGroup(basePath)) return false;
	  CHECK_ERROR(setAttributeStr("openephys:<metadata>/", basePath, "schema_id"));
	  int nMetadata = info->getEventMetadataCount();

	  timeSeries->metaDataSet.clear(); //just in case
	  for (int i = 0; i < nMetadata; i++)
	  {
		  const MetadataDescriptor* desc = info->getEventMetadataDescriptor(i);
		  String fieldName = "Field_" + String(i+1);
		  String name = desc->getName();
		  String description = desc->getDescription();
		  String identifier = desc->getIdentifier();
		  BaseDataType type = getMetadataH5Type(desc->getType(), desc->getLength()); //only string types use length, for others is always set to 1. If array types are implemented, change this
		  int length = desc->getType() == MetadataDescriptor::CHAR ? 1 : desc->getLength(); //strings are a single element of length set in the type (see above) while other elements are saved as arrays
		  String fullPath = basePath + "/" + fieldName;
		  HDF5RecordingData* dSet = createDataSet(type, 0, length, EVENT_CHUNK_SIZE, fullPath);
		  if (!dSet) return false;
		  timeSeries->metaDataSet.add(dSet);

		  CHECK_ERROR(setAttributeStr("openephys:<metadata>/", fullPath, "schema_id"));
		  CHECK_ERROR(setAttributeStr(name, fullPath, "name"));
		  CHECK_ERROR(setAttributeStr(description, fullPath, "description"));
		  CHECK_ERROR(setAttributeStr(identifier, fullPath, "identifier"));
	  }
      return true;
  }

  void NWBFile::writeEventMetadata(TimeSeries* timeSeries, const MetadataEventObject* info, const MetadataEvent* event)
  {
	  jassert(timeSeries->metaDataSet.size() == event->getMetadataValueCount());
	  jassert(info->getEventMetadataCount() == event->getMetadataValueCount());
	  int nMetadata = event->getMetadataValueCount();
	  for (int i = 0; i < nMetadata; i++)
	  {
		  BaseDataType type = getMetadataH5Type(info->getEventMetadataDescriptor(i)->getType(), info->getEventMetadataDescriptor(i)->getLength());
		  timeSeries->metaDataSet[i]->writeDataBlock(1, type, event->getMetadataValue(i)->getRawValuePointer());
	  }

  }
 
  void NWBFile::createTextDataSet(String path, String name, String text)
  {
	  ScopedPointer<HDF5RecordingData> dSet;

	  if (text.isEmpty()) text = " "; //to avoid 0-length strings, which cause errors
	  BaseDataType type = BaseDataType::STR(text.length());

	  dSet = createDataSet(type, 1, 0, path + "/" + name);
	  if (!dSet) return;
	  dSet->writeDataBlock(1, type, text.toUTF8());
  }

  void NWBFile::createBinaryDataSet(String path, String name, BaseDataType type, int length, void* data)
  {
	  ScopedPointer<HDF5RecordingData> dSet;
	  if ((length < 1) || !data) return;

	  dSet = createDataSet(type, 1, length, 1, path + "/" + name);
	  if (!dSet) return;
	  dSet->writeDataBlock(1, type, data);
  }


  //These two methods whould be easy to adapt to support array types for all base types, for now
  //length is only used for string types.
  NWBFile::BaseDataType NWBFile::getEventH5Type(EventChannel::Type type, int length)
  {
	  switch (type)
	  {
	  case EventChannel::INT8_ARRAY:
		  return BaseDataType::I8;
	  case EventChannel::UINT8_ARRAY:
		  return BaseDataType::U8;
	  case EventChannel::INT16_ARRAY:
		  return BaseDataType::I16;
	  case EventChannel::UINT16_ARRAY:
		  return BaseDataType::U16;
	  case EventChannel::INT32_ARRAY:
		  return BaseDataType::I32;
	  case EventChannel::UINT32_ARRAY:
		  return BaseDataType::U32;
	  case EventChannel::INT64_ARRAY:
		  return BaseDataType::I64;
	  case EventChannel::UINT64_ARRAY:
		  return BaseDataType::U64;
	  case EventChannel::FLOAT_ARRAY:
		  return BaseDataType::F32;
	  case EventChannel::DOUBLE_ARRAY:
		  return BaseDataType::F64;
	  case EventChannel::TEXT:
		  return BaseDataType::STR(length);
	  default:
		  return BaseDataType::I8;
	  }
  }
  NWBFile::BaseDataType NWBFile::getMetadataH5Type(MetadataDescriptor::MetadataType type, int length)
  {
	  switch (type)
	  {
	  case MetadataDescriptor::INT8:
		  return BaseDataType::I8;
	  case MetadataDescriptor::UINT8:
		  return BaseDataType::U8;
	  case MetadataDescriptor::INT16:
		  return BaseDataType::I16;
	  case MetadataDescriptor::UINT16:
		  return BaseDataType::U16;
	  case MetadataDescriptor::INT32:
		  return BaseDataType::I32;
	  case MetadataDescriptor::UINT32:
		  return BaseDataType::U32;
	  case MetadataDescriptor::INT64:
		  return BaseDataType::I64;
	  case MetadataDescriptor::UINT64:
		  return BaseDataType::U64;
	  case MetadataDescriptor::FLOAT:
		  return BaseDataType::F32;
	  case MetadataDescriptor::DOUBLE:
		  return BaseDataType::F64;
	  case MetadataDescriptor::CHAR:
		  return BaseDataType::STR(length);
	  default:
		  return BaseDataType::I8;
	  }
  }
