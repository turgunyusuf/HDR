
#include <iostream>
#include "hdr_core.h"


//#define IMAGETYPE filemapimage

int main(int argc, char ** argv)
{
	char** fileList = new char*[argc];
	fileList[0] = "exe";

	
	for(int i = 1; i<argc; i++)
	{
		fileList[i] = argv[i];
		printf("file name%d is %s",i,fileList[i]);
	}

	char *outputFileName = "hdr_image.bmp";
	printf("Image fusion is started\n");
	printf("Image numbers is %d\n", (argc-1));

	bool isColor = false;
	std::string pixelType;
	ImageResolution resolution;
	ImageResolution imResolution;
	ImageImportInfo::ICCProfile iccProfile;
	Rect2D inputUnion;


	FileNameList inputFileNameList;
	ImageImportInfo* inputInfo = NULL;
	list<ImageImportInfo*> imageInfoList;
	list<ImageImportInfo*>::iterator imageInfoIterator;

	int optind = 1;
	while (optind < (argc)) 
	{
		printf("Now reading the image%d\n",optind );
		std::string filename(fileList[optind]);
		ImageImportInfo info(filename.c_str());
		inputInfo = new ImageImportInfo(info);
		imageInfoList.push_back(inputInfo);
		inputFileNameList.push_back(filename);

		Rect2D imageROI(Point2D(inputInfo->getPosition()),
			Size2D(inputInfo->width(), inputInfo->height()));

		if (optind==1) 
		{
			inputUnion = imageROI;
			isColor = inputInfo->isColor();
			pixelType = inputInfo->getPixelType();
			resolution = ImageResolution(inputInfo->getXResolution(),
				inputInfo->getYResolution());
			iccProfile = inputInfo->getICCProfile();

		}
		else
		{
			inputUnion |= imageROI;

			if (isColor != inputInfo->isColor())
			{
				exit(1);
			}
			if (pixelType != inputInfo->getPixelType()) 
			{
				exit(1);
			}         
		}
		optind++;
	}
	printf("Image reading is done\n");
	if (resolution == ImageResolution())
	{
		imResolution = ImageResolution(DEFAULT_IMAGE_RESOLUTION,
			DEFAULT_IMAGE_RESOLUTION);
	}
	else
	{
		imResolution = resolution;
	}

	if (OutputSizeGiven)
	{
		inputUnion |= Rect2D(OutputOffsetXCmdLine,
			OutputOffsetYCmdLine,
			OutputOffsetXCmdLine + OutputWidthCmdLine,
			OutputOffsetYCmdLine + OutputHeightCmdLine);
	}
	bool rtv = 1;

	if(pixelType == "UINT8") 
	{
		printf("Image pixel type is UINT8\n");
		rtv = hdrMain<RGBValue<UInt8 > >(inputFileNameList, imageInfoList, inputUnion);
	}
	else if( pixelType == "UINT16")
	{
		printf("Image pixel type is UINT16\n");
		rtv = hdrMain<RGBValue<UInt16> >(inputFileNameList, imageInfoList,  inputUnion);
	}
	

	return 0;
}