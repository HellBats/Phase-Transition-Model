#include "include/Plots.h"
#include "include/pbPlots.h"
#include "include/supportLib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void plotLineGraph(double *xs, double *ys, int size, char *filename) {
    RGBABitmapImageReference *imageReference = CreateRGBABitmapImageReference();
    StringReference *errorMessage = CreateStringReference(L"",0);

    bool success = DrawScatterPlot(imageReference, 1080, 720, xs, size, ys, size, errorMessage);

    if (success) {
        size_t length;
        ByteArray *png = ConvertToPNG(imageReference->image);
        WriteToFile(png, filename);
    } else {
        fprintf(stderr, "Error plotting: %s\n", (char *)errorMessage->string);
    }

    FreeAllocations();
}