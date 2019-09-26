#ifndef TF_IMAGE_STATISTICS_H
#define TF_IMAGE_STATISTICS_H

#include <math.h> // powf, roundf, sqrtf
#include <stdint.h> // int32_t, uint32_t, uintptr_t

typedef int32_t tfis_s32;
typedef tfis_s32 tfis_b32;

typedef unsigned char tfis_u8;
typedef uint32_t tfis_u32;

typedef uintptr_t tfis_umm;

typedef float tfis_f32;

#define TFIS_TRUE 1
#define TFIS_FALSE 0

#define TFIS_ASSERT(Expression) if(!(Expression)) {*(int*)0 = 0;}

#define TFIS_ARRAY_COUNT(Array) (sizeof((Array)) / sizeof((Array)[0]))

#define TFIS_BOTTOM_UP_PTR(Address, Width, Height, X, Y) \
((tfis_u32*)(Address) + (X) + ((Width)*(((Height) - 1) - (Y))))

typedef struct tfis_v3
{
    tfis_f32 x, y, z;
} tfis_v3;

typedef struct tfis_v4
{
    tfis_f32 x, y, z, w;
} tfis_v4;

typedef struct tfis_random_series
{
    tfis_u32 State;
} tfis_random_series;

typedef struct tfis_counted_color
{
	tfis_u32 Color;
	int Count;
	struct tfis_counted_color* NextInHash;
} tfis_counted_color;
typedef struct tfis_color_table
{
	tfis_counted_color* Hash[1024];
} tfis_color_table;

static tfis_v3
TFIS_V3(tfis_f32 X, tfis_f32 Y, tfis_f32 Z)
{
    tfis_v3 Result;
    Result.x = X;
    Result.y = Y;
    Result.z = Z;
    
    return(Result);
}

static tfis_v4
TFIS_V4(tfis_f32 X, tfis_f32 Y, tfis_f32 Z, tfis_f32 W)
{
    tfis_v4 Result;
    Result.x = X;
    Result.y = Y;
    Result.z = Z;
    Result.w = W;
    
    return(Result);
}

// 
// NOTE(gus): Memory operations
// 

#ifndef TFIS_MALLOC
#include <stdlib.h>
#define TFIS_MALLOC(Size) malloc((Size))
#endif
static void*
TFIS_Allocate(tfis_umm Size)
{
    void* Result = TFIS_MALLOC(Size);
    
    return(Result);
}

#ifndef TFIS_FREE
#include <stdlib.h>
#define TFIS_FREE(Address) free((Address))
#endif
static void
TFIS_Free(void* Address)
{
    TFIS_FREE(Address);
}

static void
TFIS_CopyMemory(void* Dest, void* Source, tfis_umm Size)
{
    tfis_u8* DestBytePtr = (tfis_u8*)Dest;
    tfis_u8* SourceBytePtr = (tfis_u8*)Source;
    while(Size--)
    {
        *DestBytePtr++ = *SourceBytePtr++;
    }
}

#define TFIS_ZERO_ARRAY(Ptr, Count, type) TFIS_ZeroSize_((Ptr), sizeof(type)*(Count))
#define TFIS_ZERO_STRUCT(Instance, type) TFIS_ZeroSize_(&(Instance), sizeof(type))
inline void
TFIS_ZeroSize_(void* Ptr, tfis_umm Size)
{
    tfis_u8* BytePtr = (tfis_u8*)Ptr;
    while(Size--)
    {
        *BytePtr++ = 0;
    }
}

// 
// NOTE(gus): RNG
// 

static tfis_random_series
TFIS_SeedSeries(tfis_u32 Seed)
{
    tfis_random_series Result;
    
    if(Seed == 0)
    {
        Seed = 1;
    }
    
    Result.State = Seed;
    
    return(Result);
}

static tfis_u32
TFIS_GetRandom(tfis_random_series* Series)
{
    // NOTE(gus): Xorshift implementation from en.wikipedia.org/wiki/Xorshift
    
    tfis_u32 Result = Series->State;
    Result ^= (Result << 13);
	Result ^= (Result >> 17);
	Result ^= (Result << 5);
    
    Series->State = Result;
    
    return(Result);
}

static tfis_u32
TFIS_RandomU32Between(tfis_random_series* Series, tfis_u32 Min, tfis_u32 Max)
{
    // NOTE(gus): This function is inclusive of both Min and Max
    
    tfis_u32 Range = (Max - Min);
    tfis_u32 Random = TFIS_GetRandom(Series);
    tfis_u32 Result = (Min + (Random % Range));
    
    return(Result);
}

// 
// NOTE(gus): Scalar operations
// 

static tfis_f32
TFIS_Clamp(tfis_f32 Min, tfis_f32 A, tfis_f32 Max)
{
    tfis_f32 Result = A;
    if(A < Min)
    {
        Result = Min;
    }
    else if(A > Max)
    {
        Result = Max;
    }
    return(Result);
}

static tfis_f32
TFIS_Clamp01(tfis_f32 A)
{
    tfis_f32 Result = TFIS_Clamp(0.0f, A, 1.0f);
    
    return(Result);
}

static tfis_f32
TFIS_Pow(tfis_f32 X, tfis_f32 Y)
{
    tfis_f32 Result = powf(X, Y);
    
    return(Result);
}

static int
TFIS_RoundToInt(tfis_f32 A)
{
    int Result = (int)roundf(A);
    
    return(Result);
}

static tfis_u32
TFIS_RoundToU32(tfis_f32 A)
{
    tfis_u32 Result = (tfis_u32)roundf(A);
    
    return(Result);
}

static tfis_f32
TFIS_SafeRatioN(tfis_f32 Dividend, tfis_f32 Divisor, tfis_f32 N)
{
    tfis_f32 Result = N;
    if(Divisor != 0.0f)
    {
        Result = (Dividend / Divisor);
    }
    return(Result);
}

static tfis_f32
TFIS_SafeRatio0(tfis_f32 Dividend, tfis_f32 Divisor)
{
    tfis_f32 Result = TFIS_SafeRatioN(Dividend, Divisor, 0.0f);
    
    return(Result);
}

static tfis_f32
TFIS_Square(tfis_f32 A)
{
    tfis_f32 Result = (A*A);
    
    return(Result);
}

static tfis_f32
TFIS_SquareRoot(tfis_f32 A)
{
    tfis_f32 Result = sqrtf(A);
    
    return(Result);
}

static tfis_f32
TFIS_sRGBToLinear(tfis_f32 A)
{
    A = TFIS_Clamp01(A);
    
    tfis_f32 Result;
    if(A <= 0.04045f)
    {
        Result = (A / 12.92f);
    }
    else
    {
        Result = TFIS_Pow(((A + 0.055f) / 1.055f), 2.4f);
    }
    
    return(Result);
}

static tfis_f32
TFIS_LinearTosRGB(tfis_f32 A)
{
    A = TFIS_Clamp01(A);
    
    tfis_f32 Result;
    if(A <= 0.0031308f)
    {
        Result = (12.92f*A);
    }
    else
    {
        Result = ((1.055f*TFIS_Pow(A, (1.0f / 2.4f))) - 0.055f);
    }
    
    return(Result);
}

// 
// NOTE(gus): v3 operations
// 

static tfis_f32
TFIS_DotV3(tfis_v3 A, tfis_v3 B)
{
    tfis_f32 Result = ((A.x*B.x) +
                       (A.y*B.y) +
                       (A.z*B.z));
    return(Result);
}

static tfis_v3
TFIS_UnpackTriplet(tfis_u32 C)
{
    tfis_f32 Inv255 = (1.0f / 255.0f);
    tfis_v3 Result;
    Result.x = (((C >> 16) & 0xFF)*Inv255);
    Result.y = (((C >> 8) & 0xFF)*Inv255);
    Result.z = (((C >> 0) & 0xFF)*Inv255);
    
    return(Result);
}

static tfis_u32
TFIS_PackTriplet(tfis_v3 A)
{
    tfis_u32 Result = ((TFIS_RoundToU32(A.x*255.0f) << 16) |
                       (TFIS_RoundToU32(A.y*255.0f) << 8) |
                       (TFIS_RoundToU32(A.z*255.0f) << 0));
    return(Result);
}

// 
// NOTE(gus): v4 operations
// 

static tfis_v4
TFIS_sRGBToLinearV4(tfis_v4 A)
{
    tfis_v4 Result;
    Result.x = TFIS_sRGBToLinear(A.x);
    Result.y = TFIS_sRGBToLinear(A.y);
    Result.z = TFIS_sRGBToLinear(A.z);
    Result.w = TFIS_Clamp01(A.w);
    
    return(Result);
}

static tfis_v4
TFIS_UnpackRGBA(tfis_u32 C)
{
    tfis_f32 Inv255 = (1.0f / 255.0f);
    tfis_v4 Result;
    Result.x = (((C >> 0) & 0xFF)*Inv255);
    Result.y = (((C >> 8) & 0xFF)*Inv255);
    Result.z = (((C >> 16) & 0xFF)*Inv255);
    Result.w = (((C >> 24) & 0xFF)*Inv255);
    
    return(Result);
}

static tfis_u32
TFIS_PackRGBA(tfis_v4 A)
{
    tfis_u32 Result = ((TFIS_RoundToU32(A.x*255.0f) << 0) |
                       (TFIS_RoundToU32(A.y*255.0f) << 8) |
                       (TFIS_RoundToU32(A.z*255.0f) << 16) |
                       (TFIS_RoundToU32(A.w*255.0f) << 24));
    return(Result);
}

// 
// NOTE(gus): Image operations
// 

static tfis_v3
TFIS_ComputeMeanV3(void* Address, int Width, int Height)
{
    int TexelCount = 0;
    tfis_v3 Accumulator = { 0 };
    
    tfis_u32* Row = (tfis_u32*)Address;
    for(int Y = 0;
        Y < Height;
        Y++)
    {
        tfis_u32* TexelPtr = Row;
        for(int X = 0;
            X < Width;
            X++)
        {
            tfis_u32 Texel = *TexelPtr;
            tfis_v4 TexelV4 = TFIS_UnpackRGBA(Texel);
            Accumulator.x += TexelV4.x;
            Accumulator.y += TexelV4.y;
            Accumulator.z += TexelV4.z;
            
            TexelCount++;
            TexelPtr++;
        }
        
        Row += Width;
    }
    
    tfis_f32 InvTexelCount = TFIS_SafeRatio0(1.0f, (tfis_f32)TexelCount);
    tfis_v3 Result;
    Result.x = (Accumulator.x*InvTexelCount);
    Result.y = (Accumulator.y*InvTexelCount);
    Result.z = (Accumulator.z*InvTexelCount);
    
    return(Result);
}

static tfis_u32
TFIS_HashColor(tfis_u32 C)
{
    // @Refactor(gus): Better hash function please
    tfis_u32 Result = C;
    Result = (~Result + (Result << 15));
    Result = (Result ^ (Result >> 12));
    Result = (Result + (Result << 2));
    Result = (Result ^ (Result >> 4));
    Result = (Result ^ (Result >> 16));
    
    return(Result);
}

static void
TFIS_IncrementCount(tfis_color_table* Table, tfis_u32 Color)
{
    tfis_b32 EntryFound = TFIS_FALSE;
    
    tfis_u32 HashMask = (TFIS_ARRAY_COUNT(Table->Hash) - 1);
    tfis_u32 HashIndex = (TFIS_HashColor(Color) & HashMask);
    for(tfis_counted_color* Entry = Table->Hash[HashIndex];
        Entry;
        Entry = Entry->NextInHash)
    {
        if(Entry->Color == Color)
        {
            Entry->Count++;
            EntryFound = TFIS_TRUE;
            break;
        }
    }
    if(!EntryFound)
    {
        tfis_counted_color* NewEntry = (tfis_counted_color*)
            TFIS_Allocate(sizeof(tfis_counted_color));
        NewEntry->Color = Color;
        NewEntry->Count = 1;
        NewEntry->NextInHash = Table->Hash[HashIndex];
        
        Table->Hash[HashIndex] = NewEntry;
    }
}

static tfis_u32
TFIS_ComputeModeTriplet(void* Address, int Width, int Height, tfis_f32 SamplePercent)
{
    tfis_u32 Result = 0x000000;
    
    SamplePercent = TFIS_Clamp01(SamplePercent);
    
    tfis_color_table* Table = (tfis_color_table*)malloc(sizeof(tfis_color_table));
    TFIS_ZERO_STRUCT(*Table, tfis_color_table);
    
    int SampleCount = TFIS_RoundToInt((Width*Height)*SamplePercent);
    if(SampleCount)
    {
        tfis_random_series Series = TFIS_SeedSeries(2019);
        tfis_u32 SampleMask = 0xFFFCFCFC;
        for(int SampleIndex = 0;
            SampleIndex < SampleCount;
            SampleIndex++)
        {
            tfis_u32 SampleX = TFIS_RandomU32Between(&Series, 0, (tfis_u32)(Width - 1));
            tfis_u32 SampleY = TFIS_RandomU32Between(&Series, 0, (tfis_u32)(Height - 1));
            
            tfis_u32 Sample = *TFIS_BOTTOM_UP_PTR(Address, Width, Height, SampleX, SampleY);
            Sample &= SampleMask;
            
            TFIS_IncrementCount(Table, Sample);
        }
        
        int MaxCount = 0;
        tfis_u32 Mode = 0;
        for(int HashIndex = 0;
            HashIndex < TFIS_ARRAY_COUNT(Table->Hash);
            HashIndex++)
        {
            tfis_counted_color* Entry = Table->Hash[HashIndex];
            while(Entry)
            {
                if(Entry->Count > MaxCount)
                {
                    MaxCount = Entry->Count;
                    Mode = Entry->Color;
                }
                Entry = Entry->NextInHash;
            }
        }
        
        // NOTE: Repacking mode from RGBA to triplet
        Result = ((((Mode >> 0) & 0xFF) << 16) |
                  (((Mode >> 8) & 0xFF) << 8) |
                  (((Mode >> 16) & 0xFF) << 0));
        
        for(int EntryIndex = 0;
            EntryIndex < TFIS_ARRAY_COUNT(Table->Hash);
            EntryIndex++)
        {
            tfis_counted_color* Color = Table->Hash[EntryIndex];
            while(Color)
            {
                tfis_counted_color* NextInHash = Color->NextInHash;
                TFIS_Free(Color);
                Color = NextInHash;
            }
        }
    }
    
    return(Result);
}

static void
TFIS_FillBitmap(void* Address, int Width, int Height, tfis_u32 C)
{
    tfis_u32* Row = (tfis_u32*)Address;
    for(int Y = 0;
        Y < Height;
        Y++)
    {
        tfis_u32* TexelPtr = Row;
        for(int X = 0;
            X < Width;
            X++)
        {
            *TexelPtr++ = C;
        }
        
        Row += Width;
    }
}

static void
TFIS_Grayscale(void* Address, int Width, int Height)
{
    tfis_v3 Weights = TFIS_V3(0.212656f, 0.715158f, 0.072186f);
    
    tfis_u32* Row = (tfis_u32*)Address;
    for(int Y = 0;
        Y < Height;
        Y++)
    {
        tfis_u32* TexelPtr = Row;
        for(int X = 0;
            X < Width;
            X++)
        {
            tfis_v4 Texel = TFIS_UnpackRGBA(*TexelPtr);
            Texel = TFIS_sRGBToLinearV4(Texel);
            
            tfis_f32 Luminance = TFIS_DotV3(TFIS_V3(Texel.x, Texel.y, Texel.z),
                                            Weights);
            Luminance = TFIS_LinearTosRGB(Luminance);
            
            tfis_v4 GrayTexel = TFIS_V4(Luminance, Luminance, Luminance, Texel.w);
            tfis_u32 GrayTexel255 = TFIS_PackRGBA(GrayTexel);
            
            *TexelPtr++ = GrayTexel255;
        }
        
        Row += Width;
    }
}

// 
// NOTE(gus): User-facing API
// 

static void*
TFISBrightnessHistogram(void* SourceAddress, int SourceWidth, int SourceHeight,
                        int* ResultWidth, int* ResultHeight)
{
    int LuminanceBuckets[256];
    TFIS_ZERO_ARRAY(&LuminanceBuckets, TFIS_ARRAY_COUNT(LuminanceBuckets), int);
    
    tfis_umm SourceSize = (sizeof(tfis_u32)*SourceWidth*SourceHeight);
    void* GrayAddress = TFIS_Allocate(SourceSize);
    TFIS_CopyMemory(GrayAddress, SourceAddress, SourceSize);
    TFIS_Grayscale(GrayAddress, SourceWidth, SourceHeight);
    tfis_u32* Row = (tfis_u32*)GrayAddress;
    for(int Y = 0;
        Y < SourceHeight;
        Y++)
    {
        tfis_u32* TexelPtr = Row;
        for(int X = 0;
            X < SourceWidth;
            X++)
        {
            tfis_u32 Texel = *TexelPtr;
            
            tfis_u32 Luminance = (Texel & 0xFF);
            LuminanceBuckets[Luminance]++;
            
            TexelPtr++;
        }
        
        Row += SourceWidth;
    }
    TFIS_Free(GrayAddress);
    
    int HistogramWidth = TFIS_ARRAY_COUNT(LuminanceBuckets);
    if(ResultWidth) {*ResultWidth = HistogramWidth;}
    int HistogramHeight = 200;
    if(ResultHeight) {*ResultHeight = HistogramHeight;}
    
    void* Result = TFIS_Allocate(sizeof(tfis_u32)*HistogramWidth*HistogramHeight);
    TFIS_FillBitmap(Result, HistogramWidth, HistogramHeight, 0xFF000000);
    
    int MaxBucketEntries = 0;
    for(int EntryIndex = 0;
        EntryIndex < TFIS_ARRAY_COUNT(LuminanceBuckets);
        EntryIndex++)
    {
        if(LuminanceBuckets[EntryIndex] > MaxBucketEntries)
        {
            MaxBucketEntries = LuminanceBuckets[EntryIndex];
        }
    }
    tfis_f32 NormalizingC =
        (TFIS_SafeRatio0(1.0f, (tfis_f32)MaxBucketEntries)*HistogramHeight);
    
    for(int X = 0;
        X < TFIS_ARRAY_COUNT(LuminanceBuckets);
        X++)
    {
        int PixelCount = TFIS_RoundToInt(LuminanceBuckets[X]*NormalizingC);
        for(int Y = 0;
            Y < PixelCount;
            Y++)
        {
            *TFIS_BOTTOM_UP_PTR(Result, HistogramWidth, HistogramHeight, X, Y) = 0xFFFFFFFF;
        }
    }
    
    return(Result);
}

static tfis_f32
TFISCompare(void* AAddress, int AWidth, int AHeight,
            void* BAddress, int BWidth, int BHeight)
{
    tfis_f32 Result = 1.0f;
    
    if((AWidth == BWidth) &&
       (AHeight == BHeight))
    {
        tfis_f32 SquaredError = 0.0f;
        
        tfis_u32* ARow = (tfis_u32*)AAddress;
        tfis_u32* BRow = (tfis_u32*)BAddress;
        for(int Y = 0;
            Y < AHeight;
            Y++)
        {
            tfis_u32* ATexelPtr = ARow;
            tfis_u32* BTexelPtr = BRow;
            for(int X = 0;
                X < AWidth;
                X++)
            {
                tfis_u32 ATexel = *ATexelPtr;
                tfis_v4 ATexelV4 = TFIS_UnpackRGBA(ATexel);
                ATexelV4 = TFIS_sRGBToLinearV4(ATexelV4);
                
                tfis_u32 BTexel = *BTexelPtr;
                tfis_v4 BTexelV4 = TFIS_UnpackRGBA(BTexel);
                BTexelV4 = TFIS_sRGBToLinearV4(BTexelV4);
                
                tfis_f32 SquareDiffr = TFIS_Square((ATexelV4.x - BTexelV4.x)*255.0f);
                tfis_f32 SquareDiffg = TFIS_Square((ATexelV4.y - BTexelV4.y)*255.0f);
                tfis_f32 SquareDiffb = TFIS_Square((ATexelV4.z - BTexelV4.z)*255.0f);
                tfis_f32 SquareDiffa = TFIS_Square((ATexelV4.w - BTexelV4.w)*255.0f);
                
                SquaredError += SquareDiffr;
                SquaredError += SquareDiffg;
                SquaredError += SquareDiffb;
                SquaredError += SquareDiffa;
                
                ATexelPtr++;
                BTexelPtr++;
            }
            
            ARow += AWidth;
            BRow += BWidth;
        }
        
        tfis_f32 InvByteCount = (1.0f / (AWidth*AHeight*sizeof(tfis_u32)));
        tfis_f32 Inv255 = (1.0f / 255.0f);
        Result = (TFIS_SquareRoot(SquaredError*InvByteCount)*Inv255);
    }
    
    return(Result);
}

static void*
TFISHistogram(void* SourceAddress, int SourceWidth, int SourceHeight,
              int* ResultWidth, int* ResultHeight)
{
#define TFIS_BUCKET_CAPACITY 256
    int Bucketsr[TFIS_BUCKET_CAPACITY];
    int Bucketsg[TFIS_BUCKET_CAPACITY];
    int Bucketsb[TFIS_BUCKET_CAPACITY];
    TFIS_ZERO_ARRAY(&Bucketsr, TFIS_ARRAY_COUNT(Bucketsr), int);
    TFIS_ZERO_ARRAY(&Bucketsg, TFIS_ARRAY_COUNT(Bucketsg), int);
    TFIS_ZERO_ARRAY(&Bucketsb, TFIS_ARRAY_COUNT(Bucketsb), int);
    
    tfis_u32* Row = (tfis_u32*)SourceAddress;
    for(int Y = 0;
        Y < SourceHeight;
        Y++)
    {
        tfis_u32* TexelPtr = Row;
        for(int X = 0;
            X < SourceWidth;
            X++)
        {
            tfis_u32 Texel = *TexelPtr;
            
            tfis_u32 Texelr = ((Texel >> 0) & 0xFF);
            tfis_u32 Texelg = ((Texel >> 8) & 0xFF);
            tfis_u32 Texelb = ((Texel >> 16) & 0xFF);
            
            Bucketsr[Texelr]++;
            Bucketsg[Texelg]++;
            Bucketsb[Texelb]++;
            
            TexelPtr++;
        }
        
        Row += SourceWidth;
    }
    
    int HistogramWidth = TFIS_BUCKET_CAPACITY;
    if(ResultWidth) {*ResultWidth = HistogramWidth;}
    int HistogramHeight = 200;
    if(ResultHeight) {*ResultHeight = HistogramHeight;}
    
    void* Result = TFIS_Allocate(sizeof(tfis_u32)*HistogramWidth*HistogramHeight);
    TFIS_FillBitmap(Result, HistogramWidth, HistogramHeight, 0xFF000000);
    
    int MaxBucketEntries = 0;
    for(int EntryIndex = 0;
        EntryIndex < TFIS_BUCKET_CAPACITY;
        EntryIndex++)
    {
        if(Bucketsr[EntryIndex] > MaxBucketEntries)
        {
            MaxBucketEntries = Bucketsr[EntryIndex];
        }
        if(Bucketsg[EntryIndex] > MaxBucketEntries)
        {
            MaxBucketEntries = Bucketsg[EntryIndex];
        }
        if(Bucketsb[EntryIndex] > MaxBucketEntries)
        {
            MaxBucketEntries = Bucketsb[EntryIndex];
        }
    }
    tfis_f32 NormalizingC =
        (TFIS_SafeRatio0(1.0f, (tfis_f32)MaxBucketEntries)*HistogramHeight);
    
    for(int X = 0;
        X < TFIS_BUCKET_CAPACITY;
        X++)
    {
        int PixelCountr = TFIS_RoundToInt(Bucketsr[X]*NormalizingC);
        int PixelCountg = TFIS_RoundToInt(Bucketsg[X]*NormalizingC);
        int PixelCountb = TFIS_RoundToInt(Bucketsb[X]*NormalizingC);
        
        for(int Y = 0;
            Y < PixelCountr;
            Y++)
        {
            *TFIS_BOTTOM_UP_PTR(Result, HistogramWidth, HistogramHeight, X, Y) |= 0xFF0000FF;
        }
        for(int Y = 0;
            Y < PixelCountg;
            Y++)
        {
            *TFIS_BOTTOM_UP_PTR(Result, HistogramWidth, HistogramHeight, X, Y) |= 0xFF00FF00;
        }
        for(int Y = 0;
            Y < PixelCountb;
            Y++)
        {
            *TFIS_BOTTOM_UP_PTR(Result, HistogramWidth, HistogramHeight, X, Y) |= 0xFFFF0000;
        }
    }
#undef TFIS_BUCKET_CAPACITY
    
    return(Result);
}

static void
TFISMeanFloat(void* SourceAddress, int SourceWidth, int SourceHeight,
              tfis_f32* R, tfis_f32* G, tfis_f32* B)
{
    tfis_v3 MeanV3 = TFIS_ComputeMeanV3(SourceAddress, SourceWidth, SourceHeight);
    *R = MeanV3.x;
    *G = MeanV3.y;
    *B = MeanV3.z;
}

static void
TFISMeanTriplet(void* SourceAddress, int SourceWidth, int SourceHeight,
                tfis_u32* Triplet)
{
    tfis_v3 MeanV3 = TFIS_ComputeMeanV3(SourceAddress, SourceWidth, SourceHeight);
    tfis_u32 Result = TFIS_PackTriplet(MeanV3);
    *Triplet = Result;
}

static void
TFISModeFloat(void* SourceAddress, int SourceWidth, int SourceHeight,
              tfis_f32 SamplePercent,
              tfis_f32* R, tfis_f32* G, tfis_f32* B)
{
    tfis_u32 Mode = TFIS_ComputeModeTriplet(SourceAddress, SourceWidth, SourceHeight,
                                            SamplePercent);
    tfis_v3 ModeV3 = TFIS_UnpackTriplet(Mode);
    *R = ModeV3.x;
    *G = ModeV3.y;
    *B = ModeV3.z;
}

static void
TFISModeTriplet(void* SourceAddress, int SourceWidth, int SourceHeight,
                tfis_f32 SamplePercent,
                tfis_u32* Triplet)
{
    tfis_u32 Mode = TFIS_ComputeModeTriplet(SourceAddress, SourceWidth, SourceHeight,
                                            SamplePercent);
    *Triplet = Mode;
}

static void
TFISFree(void* Address)
{
    TFIS_Free(Address);
}

#endif

/*
MIT License

Copyright (c) 2019 Gustavo Velasquez (@twelvefifteen on GitHub)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this softwareand associated documentation files(the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and /or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions :

The above copyright noticeand this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
