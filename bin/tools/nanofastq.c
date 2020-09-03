#pragma warning(disable : 4996)

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <zlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <assert.h>
#include <math.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define PHRED_0_VALUE	33

int main(int argc, const char **argv)
{
    // Construct usageText
#define USAGE_TEXT_MAX_LENGTH	1024
	char usageText[USAGE_TEXT_MAX_LENGTH];
	{
		memset(usageText, ' ', USAGE_TEXT_MAX_LENGTH);

#define USAGE_TEXT_NUM_LINE		6
		char usageTextLine0[] = "Usage: input | %s <options> | output\n";
		char usageTextLine1[] = "  options: -q minimum average read quality\n";
		char usageTextLine2[] = "           -l minimum read length\n";
		char usageTextLine3[] = "           -h headcrop\n";
		char usageTextLine4[] = "           -t tailcrop\n";
		char usageTextLine5[] = "           -r read ID prefix (for reassigning read_id)";
		char *usageTextArray[USAGE_TEXT_NUM_LINE];
		usageTextArray[0] = usageTextLine0;
		usageTextArray[1] = usageTextLine1;
		usageTextArray[2] = usageTextLine2;
		usageTextArray[3] = usageTextLine3;
		usageTextArray[4] = usageTextLine4;
		usageTextArray[5] = usageTextLine5;

		int usageTextOffset = 0;
		for (int i = 0; i < USAGE_TEXT_NUM_LINE; ++i)
		{
			strcpy(usageText + usageTextOffset, usageTextArray[i]);
			usageTextOffset += (int)strlen(usageTextArray[i]);
			usageText[usageTextOffset] = ' ';
		}
		usageText[usageTextOffset] = '\0';
	}

	float minQuality = 0;
	int minLength = 0;
	int headcrop = 0;
	int tailcrop = 0;
	int reassignReadID = 0;
	char *readID_prefix = NULL;

	int argumentCount = 0;
	while (argumentCount + 1 < argc)
	{
		if (strcmp(argv[argumentCount + 1], "-q") == 0)
		{
			if (argumentCount + 1 >= argc)
			{
				fprintf(stderr, usageText, argv[0]);
				exit(EXIT_FAILURE);
			}
			minQuality = strtof(argv[argumentCount + 2],NULL);
			argumentCount += 2;
			continue;
		}
		if (strcmp(argv[argumentCount + 1], "-l") == 0)
		{
			if (argumentCount + 1 >= argc)
			{
				fprintf(stderr, usageText, argv[0]);
				exit(EXIT_FAILURE);
			}
			minLength = atoi(argv[argumentCount + 2]);
			argumentCount += 2;
			continue;
		}
		if (strcmp(argv[argumentCount + 1], "-h") == 0)
		{
			if (argumentCount + 1 >= argc)
			{
				fprintf(stderr, usageText, argv[0]);
				exit(EXIT_FAILURE);
			}
			headcrop = atoi(argv[argumentCount + 2]);
			argumentCount += 2;
			continue;
		}
		if (strcmp(argv[argumentCount + 1], "-t") == 0)
		{
			if (argumentCount + 1 >= argc)
			{
				fprintf(stderr, usageText, argv[0]);
				exit(EXIT_FAILURE);
			}
			tailcrop = atoi(argv[argumentCount + 2]);
			argumentCount += 2;
			continue;
		}
		if (strcmp(argv[argumentCount + 1], "-r") == 0)
		{
			if (argumentCount + 1 >= argc)
			{
				fprintf(stderr, usageText, argv[0]);
				exit(EXIT_FAILURE);
			}
			reassignReadID = 1;
			readID_prefix = (char*)argv[argumentCount + 2];
			argumentCount += 2;
			continue;
		}
		break;
	}

	if (minQuality < 0 || minLength < 0 || headcrop < 0 || tailcrop < 0)
	{
		fprintf(stderr, usageText, argv[0]);
		exit(EXIT_FAILURE);
	}
	if (minLength == 0)
	{
		minLength = 1;
	}


	gzFile inputFile;
	kseq_t *seq;

	inputFile = gzdopen(fileno(stdin), "r");
	if (inputFile == 0) {
		exit(EXIT_FAILURE);
	}

	seq = kseq_init(inputFile);

	double PhredScoreToErrorProbability[128];
	for (int i = 0; i < 128; ++i)
		PhredScoreToErrorProbability[i] = pow(10., -i / 10.);

	char *readID;
	char *serial;
	if (reassignReadID == 1)
	{
                int serialPos = strlen(readID_prefix);
		readID = (char*)calloc(serialPos + 100i, 1);
		serial = readID + serialPos;
		strcpy(readID, readID_prefix);
	}

        int numRead = 0;

	while (kseq_read(seq) >= 0) {

                numRead++;
		//char *readName;
		if (reassignReadID == 0)
		{
			readID = seq->name.s;
		}
		else
		{
			sprintf(serial, "%d", numRead);
		}

		int passedFilter = 1;

		double totalErrorProbability = 0;
		double averagePhredScore = 0;
		if (seq->is_fastq)
		{
			for (int i = 0; i < seq->seq.l; ++i)
			{
				totalErrorProbability += PhredScoreToErrorProbability[seq->qual.s[i] - PHRED_0_VALUE];
			}
			averagePhredScore = -10 * log10(totalErrorProbability / seq->seq.l);
		}

		double totalErrorProbability_afterCrop = totalErrorProbability;
		double averagePhredScore_afterCrop = 0;

		int startPos = headcrop;
		int endPos = seq->seq.l - tailcrop;

		if (endPos - startPos < minLength)
		{
			passedFilter = 0;
		}
		else
		{
			if (seq->is_fastq)
			{
				for (int i = 0; i < startPos; ++i)
				{
					totalErrorProbability_afterCrop -= PhredScoreToErrorProbability[seq->qual.s[i] - PHRED_0_VALUE];
				}
				for (int i = endPos; i < seq->seq.l; ++i)
				{
					totalErrorProbability_afterCrop -= PhredScoreToErrorProbability[seq->qual.s[i] - PHRED_0_VALUE];
				}
				averagePhredScore_afterCrop = -10 * log10(totalErrorProbability_afterCrop / (seq->seq.l - headcrop - tailcrop));
			}

			if (seq->is_fastq && averagePhredScore_afterCrop < minQuality)
			{
				passedFilter = 0;
			}
			else
			{
				putchar(seq->is_fastq ? '@' : '>');
				fputs(readID, stdout);
				if (seq->comment.l) {
					putchar(' '); puts(seq->comment.s);
				}
				else
				{
					putchar('\n');
				}
				fwrite(seq->seq.s + startPos, 1, endPos - startPos, stdout); putchar('\n');
				if (seq->is_fastq) {
					puts("+");
					fwrite(seq->qual.s + startPos, 1, endPos - startPos, stdout); putchar('\n');
				}

			}
		}

		// Output read info for statistics generation
		fputs(readID, stderr);
		fprintf(stderr, "\t%zu", seq->seq.l);
		fprintf(stderr, "\t%.2f", averagePhredScore);
		if (seq->seq.l - headcrop - tailcrop > 0)
			fprintf(stderr, "\t%zu", seq->seq.l - headcrop - tailcrop);
		else
			fprintf(stderr, "\t%d", 0);
		fprintf(stderr, "\t%.2f", averagePhredScore_afterCrop);
		fprintf(stderr, "\t%d", passedFilter);
		fputc('\n', stderr);

	}

	kseq_destroy(seq);
	gzclose(inputFile);

	return(EXIT_SUCCESS);

}
