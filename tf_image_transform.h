// ============================================================================
//  
//    tf_image_transform v0.1a - Image transforms by @twelvefifteen on GitHub
// 
// USAGE 
// 
//    The user-facing functions can be found in the bottommost code section of
//    this file.
//
//    void* ResizedAddress = TFITResize(ImageAddress, 256, 512,
//                                      128, 256,
//                                      TFITFilterType_Bilinear);
//    // This call resizes 256x512 image at ImageAddress to 128x256 and returns
//    // the raster.
//    // Resizing here uses bilinear filtering. An explanation of the available
//    // available filters can be found in the FILTERING DETAILS section below.
//    TFITFree(ResizedAddress);
// 
//    void* ClockwiseAddress = TFITRotate(ImageAddress, 256, 512,
//                                        TFITRotationDir_Clockwise);
//    // This call rotates the 256x512 image at ImageAddress clockwise.
//    // The dimensions of the rotated image will always be the reverse of the
//    // old dimensions (512x256 in this example).
//    // Counterclockwise rotation is also supported.
//    TFITFree(ClockwiseAddress);
// 
// INTERNALS
//    
//    All raster data passed to this library must contain four components and
//    be RGBA ordered.
// 
//    Allocations are made using the TFIT_MALLOC macro, which defaults to
//    malloc(). You can, however, #define this macro to use your own allocation
//    scheme. The same goes for TFIT_FREE, which defaults to free().
// 
// FILTERING DETAILS
// 
//    There are three filtering options for the resize function. Smoothness of
//    the resized images improves from nearest neighbor to bilinear to bicubic
//    filtering. However, the cost of computation also increases in this order.
// 
//    This library currently only supports Catmull-Rom splines for bicubic
//    filtering.
//  
// LICENSE
// 
//    The license for this library can be found at the bottom of this file.
// 
// ============================================================================

#ifndef TF_IMAGE_TRANSFORM_H
#define TF_IMAGE_TRANSFORM_H

#include <math.h> // powf, roundf
#include <stdint.h>

typedef unsigned char tfit_u8;
typedef uint32_t tfit_u32;

typedef uintptr_t tfit_umm;

typedef float tfit_f32;

typedef enum tfit_filter_type
{
    TFITFilterType_NearestNeighbor,
    TFITFilterType_Bilinear,
    TFITFilterType_Bicubic,
    
    TFITFilterType_Count,
} tfit_filter_type;

typedef enum tfit_rotation_direction
{
    TFITRotationDir_Clockwise,
    TFITRotationDir_Counterclockwise,
} tfit_rotation_direction;

typedef struct tfit_v4
{
    tfit_f32 x, y, z, w;
} tfit_v4;

#define TFIT_GET_BITMAP_PTR(Bitmap, MinX, MinY) ((tfit_u8*)(&(Bitmap))->Address + ((MinX)*4) + ((MinY)*(Bitmap).Pitch))
typedef struct tfit_bitmap
{
    void* Address;
    int Width;
    int Height;
    int Pitch;
} tfit_bitmap;

#define TFIT_SAMPLE_BITMAP(Name) tfit_u32 Name(tfit_bitmap Bitmap, tfit_f32 U, tfit_f32 V)
typedef TFIT_SAMPLE_BITMAP(tfit_sample_bitmap);

static tfit_bitmap
TFIT_MakeBitmap(void* Address, int Width, int Height)
{
    tfit_bitmap Result;
    Result.Address = Address;
    Result.Width = Width;
    Result.Height = Height;
    Result.Pitch = (int)(Width*sizeof(tfit_u32));
    
    return(Result);
}

// 
// Memory management
// 

#ifndef TFIT_MALLOC
#include <stdio.h>
#define TFIT_MALLOC(Size) malloc((Size))
#endif
static void*
TFIT_Allocate(tfit_umm Size)
{
    void* Result = TFIT_MALLOC(Size);
    
    return(Result);
}

#ifndef TFIT_FREE
#include <stdio.h>
#define TFIT_FREE(Address) free((Address))
#endif
static void
TFIT_Free(void* Address)
{
    TFIT_FREE(Address);
}

// 
// Scalar operations
// 

static tfit_f32 TFIT_Cube(tfit_f32 A);
static tfit_f32 TFIT_Square(tfit_f32 A);
static tfit_f32
TFIT_CatmullRomInterpolate(tfit_f32 FPrime0, tfit_f32 F0, tfit_f32 F1, tfit_f32 FPrime1, tfit_f32 t)
{
    tfit_f32 CubicC = ((-0.5f*FPrime0) + (1.5f*F0) - (1.5f*F1) + (0.5f*FPrime1));
    tfit_f32 QuadraticC = (FPrime0 + (-2.5f*F0) + (2.0f*F1) - (0.5f*FPrime1));
    tfit_f32 LinearC = ((-0.5f*FPrime0) + (0.5f*F1));
    tfit_f32 Constant = F0;
    
    tfit_f32 Result =
        ((CubicC*TFIT_Cube(t)) + (QuadraticC*TFIT_Square(t)) + (LinearC*t) + Constant);
    
    return(Result);
}

static tfit_f32
TFIT_Clamp(tfit_f32 Min, tfit_f32 A, tfit_f32 Max)
{
    tfit_f32 Result = A;
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

static int
TFIT_Clampi(int Min, int A, int Max)
{
    int Result = (int)TFIT_Clamp((tfit_f32)Min, (tfit_f32)A, (tfit_f32)Max);
    
    return(Result);
}

static tfit_f32
TFIT_Clamp01(tfit_f32 A)
{
    tfit_f32 Result = TFIT_Clamp(0.0f, A, 1.0f);
    
    return(Result);
}

static tfit_f32
TFIT_Cube(tfit_f32 A)
{
    tfit_f32 Result = (A*A*A);
    
    return(Result);
}

static tfit_f32
TFIT_Lerp(tfit_f32 A, tfit_f32 t, tfit_f32 B)
{
    tfit_f32 Result = (((1.0f - t)*A) + (t*B));
    
    return(Result);
}

static tfit_f32
TFIT_Pow(tfit_f32 X, tfit_f32 Y)
{
    tfit_f32 Result = powf(X, Y);
    
    return(Result);
}

static tfit_u32
TFIT_RoundToU32(tfit_f32 A)
{
    tfit_u32 Result = (tfit_u32)roundf(A);
    
    return(Result);
}

static tfit_f32
TFIT_SafeRatioN(tfit_f32 Dividend, tfit_f32 Divisor, tfit_f32 N)
{
    tfit_f32 Result = N;
    if(Divisor != 0.0f)
    {
        Result = (Dividend / Divisor);
    }
    return(Result);
}

static tfit_f32
TFIT_SafeRatio0(tfit_f32 Dividend, tfit_f32 Divisor)
{
    tfit_f32 Result = TFIT_SafeRatioN(Dividend, Divisor, 0.0f);
    
    return(Result);
}

static tfit_f32
TFIT_Square(tfit_f32 A)
{
    tfit_f32 Result = (A*A);
    
    return(Result);
}

static tfit_f32
TFIT_sRGBToLinear(tfit_f32 A)
{
    A = TFIT_Clamp01(A);
    
    tfit_f32 Result;
    if(A <= 0.04045f)
    {
        Result = (A / 12.92f);
    }
    else
    {
        Result = TFIT_Pow(((A + 0.055f) / 1.055f), 2.4f);
    }
    
    return(Result);
}

static tfit_f32
TFIT_LinearTosRGB(tfit_f32 A)
{
    A = TFIT_Clamp01(A);
    
    tfit_f32 Result;
    if(A <= 0.0031308f)
    {
        Result = (12.92f*A);
    }
    else
    {
        Result = ((1.055f*TFIT_Pow(A, (1.0f / 2.4f))) - 0.055f);
    }
    
    return(Result);
}

// 
// v4 operations
// 

static tfit_v4
TFIT_CatmullRomInterpolateV4(tfit_v4 FPrime0, tfit_v4 F0, tfit_v4 F1, tfit_v4 FPrime1, tfit_f32 t)
{
    tfit_v4 Result;
    Result.x = TFIT_CatmullRomInterpolate(FPrime0.x, F0.x, F1.x, FPrime1.x, t);
    Result.y = TFIT_CatmullRomInterpolate(FPrime0.y, F0.y, F1.y, FPrime1.y, t);
    Result.z = TFIT_CatmullRomInterpolate(FPrime0.z, F0.z, F1.z, FPrime1.z, t);
    Result.w = TFIT_CatmullRomInterpolate(FPrime0.w, F0.w, F1.w, FPrime1.w, t);
    
    return(Result);
}

static tfit_v4
TFIT_UnpackRGBA(tfit_u32 C)
{
    tfit_f32 Inv255 = (1.0f / 255.0f);
    tfit_v4 Result;
    Result.x = (((C >> 0) & 0xFF)*Inv255);
    Result.y = (((C >> 8) & 0xFF)*Inv255);
    Result.z = (((C >> 16) & 0xFF)*Inv255);
    Result.w = (((C >> 24) & 0xFF)*Inv255);
    
    return(Result);
}

static tfit_v4
TFIT_sRGBToLinearV4(tfit_v4 A)
{
    tfit_v4 Result;
    Result.x = TFIT_sRGBToLinear(A.x);
    Result.y = TFIT_sRGBToLinear(A.y);
    Result.z = TFIT_sRGBToLinear(A.z);
    Result.w = TFIT_Clamp01(A.w);
    
    return(Result);
}

static tfit_v4
TFIT_LinearTosRGBV4(tfit_v4 A)
{
    tfit_v4 Result;
    Result.x = TFIT_LinearTosRGB(A.x);
    Result.y = TFIT_LinearTosRGB(A.y);
    Result.z = TFIT_LinearTosRGB(A.z);
    Result.w = TFIT_Clamp01(A.w);
    
    return(Result);
}

static tfit_v4
TFIT_LerpV4(tfit_v4 A, tfit_f32 t, tfit_v4 B)
{
    tfit_v4 Result;
    Result.x = TFIT_Lerp(A.x, t, B.x);
    Result.y = TFIT_Lerp(A.y, t, B.y);
    Result.z = TFIT_Lerp(A.z, t, B.z);
    Result.w = TFIT_Lerp(A.w, t, B.w);
    
    return(Result);
}

static tfit_u32
TFIT_PackRGBA(tfit_v4 A)
{
    tfit_u32 Result = ((TFIT_RoundToU32(A.x*255.0f) << 0) |
                       (TFIT_RoundToU32(A.y*255.0f) << 8) |
                       (TFIT_RoundToU32(A.z*255.0f) << 16) |
                       (TFIT_RoundToU32(A.w*255.0f) << 24));
    return(Result);
}

// 
// Samplers
// 

static TFIT_SAMPLE_BITMAP(TFIT_SampleBitmapNearestNeighbor)
{
    tfit_u32 SampleX = TFIT_RoundToU32(U*((tfit_f32)Bitmap.Width - 1.0f));
    tfit_u32 SampleY = TFIT_RoundToU32(V*((tfit_f32)Bitmap.Height - 1.0f));
    
    tfit_u32 Result = *(tfit_u32*)TFIT_GET_BITMAP_PTR(Bitmap, SampleX, SampleY);
    
    return(Result);
}

static TFIT_SAMPLE_BITMAP(TFIT_SampleBitmapBilinear)
{
    tfit_f32 tSampleX = (U*(Bitmap.Width - 1.0f));
    tfit_f32 tSampleY = (V*(Bitmap.Height - 1.0f));
    
    int SampleX = (int)tSampleX;
    int SampleY = (int)tSampleY;
    
    tfit_f32 fSampleX = (tSampleX - SampleX);
    tfit_f32 fSampleY = (tSampleY - SampleY);
    
    int TexelAX = SampleX;
    int TexelAY = SampleY;
    
    int TexelBX = (SampleX + 1);
    int TexelBY = SampleY;
    
    int TexelCX = SampleX;
    int TexelCY = (SampleY + 1);
    
    int TexelDX = (SampleX + 1);
    int TexelDY = (SampleY + 1);
    
    // Repeating edge texels out to infinity
    if(SampleX == (Bitmap.Width - 1))
    {
        TexelBX = SampleX;
        TexelDX = SampleX;
    }
    if(SampleY == (Bitmap.Height - 1))
    {
        TexelCY = SampleY;
        TexelDY = SampleY;
    }
    
    // TexelA
    
    tfit_u32 TexelA = *(tfit_u32*)TFIT_GET_BITMAP_PTR(Bitmap, TexelAX, TexelAY);
    tfit_v4 vTexelA = TFIT_UnpackRGBA(TexelA);
    vTexelA = TFIT_sRGBToLinearV4(vTexelA);
    
    // TexelB
    
    tfit_u32 TexelB = *(tfit_u32*)TFIT_GET_BITMAP_PTR(Bitmap, TexelBX, TexelBY);
    tfit_v4 vTexelB = TFIT_UnpackRGBA(TexelB);
    vTexelB = TFIT_sRGBToLinearV4(vTexelB);
    
    // TexelC
    
    tfit_u32 TexelC = *(tfit_u32*)TFIT_GET_BITMAP_PTR(Bitmap, TexelCX, TexelCY);
    tfit_v4 vTexelC = TFIT_UnpackRGBA(TexelC);
    vTexelC = TFIT_sRGBToLinearV4(vTexelC);
    
    // TexelD
    
    tfit_u32 TexelD = *(tfit_u32*)TFIT_GET_BITMAP_PTR(Bitmap, TexelDX, TexelDY);
    tfit_v4 vTexelD = TFIT_UnpackRGBA(TexelD);
    vTexelD = TFIT_sRGBToLinearV4(vTexelD);
    
    // 
    
    tfit_v4 vResult = TFIT_LerpV4(TFIT_LerpV4(vTexelA, fSampleX, vTexelB),
                                  fSampleY,
                                  TFIT_LerpV4(vTexelC, fSampleX, vTexelD));
    vResult = TFIT_LinearTosRGBV4(vResult);
    tfit_u32 Result = TFIT_PackRGBA(vResult);
    
    return(Result);
}

static TFIT_SAMPLE_BITMAP(TFIT_SampleBitmapBicubic)
{
    tfit_f32 tSampleX = (U*(Bitmap.Width - 1.0f));
    tfit_f32 tSampleY = (V*(Bitmap.Height - 1.0f));
    
    int SampleX = (int)tSampleX;
    int SampleY = (int)tSampleY;
    
    tfit_f32 fSampleX = (tSampleX - SampleX);
    tfit_f32 fSampleY = (tSampleY - SampleY);
    
    int FPrime0X = (SampleX - 1);
    int F0X = SampleX;
    int F1X = (SampleX + 1);
    int FPrime1X = (SampleX + 2);
    
    FPrime0X = TFIT_Clampi(0, FPrime0X, (Bitmap.Width - 1));
    F0X = TFIT_Clampi(0, F0X, (Bitmap.Width - 1));
    F1X = TFIT_Clampi(0, F1X, (Bitmap.Width - 1));
    FPrime1X = TFIT_Clampi(0, FPrime1X, (Bitmap.Width - 1));
    
    int TexelYA = (SampleY - 1);
    int TexelYB = SampleY;
    int TexelYC = (SampleY + 1);
    int TexelYD = (SampleY + 2);
    
    TexelYA = TFIT_Clampi(0, TexelYA, (Bitmap.Height - 1));
    TexelYB = TFIT_Clampi(0, TexelYB, (Bitmap.Height - 1));
    TexelYC = TFIT_Clampi(0, TexelYC, (Bitmap.Height - 1));
    TexelYD = TFIT_Clampi(0, TexelYD, (Bitmap.Height - 1));
    
    // A
    
    tfit_u32 FPrime0A = *(tfit_u32*)TFIT_GET_BITMAP_PTR(Bitmap, FPrime0X, TexelYA);
    tfit_u32 F0A = *(tfit_u32*)TFIT_GET_BITMAP_PTR(Bitmap, F0X, TexelYA);
    tfit_u32 F1A = *(tfit_u32*)TFIT_GET_BITMAP_PTR(Bitmap, F1X, TexelYA);
    tfit_u32 FPrime1A = *(tfit_u32*)TFIT_GET_BITMAP_PTR(Bitmap, FPrime1X, TexelYA);
    
    tfit_v4 vFPrime0A = TFIT_UnpackRGBA(FPrime0A);
    tfit_v4 vF0A = TFIT_UnpackRGBA(F0A);
    tfit_v4 vF1A = TFIT_UnpackRGBA(F1A);
    tfit_v4 vFPrime1A = TFIT_UnpackRGBA(FPrime1A);
    
    vFPrime0A = TFIT_sRGBToLinearV4(vFPrime0A);
    vF0A = TFIT_sRGBToLinearV4(vF0A);
    vF1A = TFIT_sRGBToLinearV4(vF1A);
    vFPrime1A = TFIT_sRGBToLinearV4(vFPrime1A);
    
    tfit_v4 vBlendedA = TFIT_CatmullRomInterpolateV4(vFPrime0A, vF0A, vF1A, vFPrime1A,
                                                     fSampleX);
    
    // B
    
    tfit_u32 FPrime0B = *(tfit_u32*)TFIT_GET_BITMAP_PTR(Bitmap, FPrime0X, TexelYB);
    tfit_u32 F0B = *(tfit_u32*)TFIT_GET_BITMAP_PTR(Bitmap, F0X, TexelYB);
    tfit_u32 F1B = *(tfit_u32*)TFIT_GET_BITMAP_PTR(Bitmap, F1X, TexelYB);
    tfit_u32 FPrime1B = *(tfit_u32*)TFIT_GET_BITMAP_PTR(Bitmap, FPrime1X, TexelYB);
    
    tfit_v4 vFPrime0B = TFIT_UnpackRGBA(FPrime0B);
    tfit_v4 vF0B = TFIT_UnpackRGBA(F0B);
    tfit_v4 vF1B = TFIT_UnpackRGBA(F1B);
    tfit_v4 vFPrime1B = TFIT_UnpackRGBA(FPrime1B);
    
    vFPrime0B = TFIT_sRGBToLinearV4(vFPrime0B);
    vF0B = TFIT_sRGBToLinearV4(vF0B);
    vF1B = TFIT_sRGBToLinearV4(vF1B);
    vFPrime1B = TFIT_sRGBToLinearV4(vFPrime1B);
    
    tfit_v4 vBlendedB = TFIT_CatmullRomInterpolateV4(vFPrime0B, vF0B, vF1B, vFPrime1B,
                                                     fSampleX);
    
    // C
    
    tfit_u32 FPrime0C = *(tfit_u32*)TFIT_GET_BITMAP_PTR(Bitmap, FPrime0X, TexelYC);
    tfit_u32 F0C = *(tfit_u32*)TFIT_GET_BITMAP_PTR(Bitmap, F0X, TexelYC);
    tfit_u32 F1C = *(tfit_u32*)TFIT_GET_BITMAP_PTR(Bitmap, F1X, TexelYC);
    tfit_u32 FPrime1C = *(tfit_u32*)TFIT_GET_BITMAP_PTR(Bitmap, FPrime1X, TexelYC);
    
    tfit_v4 vFPrime0C = TFIT_UnpackRGBA(FPrime0C);
    tfit_v4 vF0C = TFIT_UnpackRGBA(F0C);
    tfit_v4 vF1C = TFIT_UnpackRGBA(F1C);
    tfit_v4 vFPrime1C = TFIT_UnpackRGBA(FPrime1C);
    
    vFPrime0C = TFIT_sRGBToLinearV4(vFPrime0C);
    vF0C = TFIT_sRGBToLinearV4(vF0C);
    vF1C = TFIT_sRGBToLinearV4(vF1C);
    vFPrime1C = TFIT_sRGBToLinearV4(vFPrime1C);
    
    tfit_v4 vBlendedC = TFIT_CatmullRomInterpolateV4(vFPrime0C, vF0C, vF1C, vFPrime1C,
                                                     fSampleX);
    
    // D
    
    tfit_u32 FPrime0D = *(tfit_u32*)TFIT_GET_BITMAP_PTR(Bitmap, FPrime0X, TexelYD);
    tfit_u32 F0D = *(tfit_u32*)TFIT_GET_BITMAP_PTR(Bitmap, F0X, TexelYD);
    tfit_u32 F1D = *(tfit_u32*)TFIT_GET_BITMAP_PTR(Bitmap, F1X, TexelYD);
    tfit_u32 FPrime1D = *(tfit_u32*)TFIT_GET_BITMAP_PTR(Bitmap, FPrime1X, TexelYD);
    
    tfit_v4 vFPrime0D = TFIT_UnpackRGBA(FPrime0D);
    tfit_v4 vF0D = TFIT_UnpackRGBA(F0D);
    tfit_v4 vF1D = TFIT_UnpackRGBA(F1D);
    tfit_v4 vFPrime1D = TFIT_UnpackRGBA(FPrime1D);
    
    vFPrime0D = TFIT_sRGBToLinearV4(vFPrime0D);
    vF0D = TFIT_sRGBToLinearV4(vF0D);
    vF1D = TFIT_sRGBToLinearV4(vF1D);
    vFPrime1D = TFIT_sRGBToLinearV4(vFPrime1D);
    
    tfit_v4 vBlendedD = TFIT_CatmullRomInterpolateV4(vFPrime0D, vF0D, vF1D, vFPrime1D,
                                                     fSampleX);
    
    // 
    
    tfit_v4 vResult = TFIT_CatmullRomInterpolateV4(vBlendedA, vBlendedB, vBlendedC, vBlendedD,
                                                   fSampleY);
    vResult = TFIT_LinearTosRGBV4(vResult);
    tfit_u32 Result = TFIT_PackRGBA(vResult);
    
    return(Result);
}

// 
// User-facing API
// 

static void*
TFITResize(void* SourceAddress,
           int SourceWidth, int SourceHeight,
           int ResizedWidth, int ResizedHeight,
           tfit_filter_type FilterType)
{
    int MinX = 0;
    int MinY = 0;
    int MaxX = ResizedWidth;
    int MaxY = ResizedHeight;
    
    tfit_sample_bitmap* SampleBitmap = 0;
    switch(FilterType)
    {
        case TFITFilterType_NearestNeighbor:
        {
            SampleBitmap = TFIT_SampleBitmapNearestNeighbor;
        } break;
        
        case TFITFilterType_Bilinear:
        {
            SampleBitmap = TFIT_SampleBitmapBilinear;
        } break;
        
        case TFITFilterType_Bicubic:
        {
            SampleBitmap = TFIT_SampleBitmapBicubic;
        } break;
        
        case TFITFilterType_Count:
        {
            
        } break;
    }
    
    tfit_bitmap SourceBitmap = TFIT_MakeBitmap(SourceAddress, SourceWidth, SourceHeight);
    void* Result = TFIT_Allocate(sizeof(tfit_u32)*ResizedWidth*ResizedHeight);
    if(Result)
    {
        tfit_u32* DestRow = (tfit_u32*)Result;
        for(int Y = MinY;
            Y < MaxY;
            Y++)
        {
            tfit_u32* DestTexelPtr = DestRow;
            for(int X = MinX;
                X < MaxX;
                X++)
            {
                tfit_f32 U = TFIT_SafeRatio0((tfit_f32)X, ((tfit_f32)MaxX - 1.0f));
                tfit_f32 V = TFIT_SafeRatio0((tfit_f32)Y, ((tfit_f32)MaxY - 1.0f));
                
                tfit_u32 Sample = SampleBitmap(SourceBitmap, U, V);
                *DestTexelPtr++ = Sample;
            }
            
            DestRow += ResizedWidth;
        }
    }
    
    return(Result);
}

static void*
TFITRotate(void* SourceAddress,
           int SourceWidth, int SourceHeight,
           tfit_rotation_direction Direction)
{
    int MinX = 0;
    int MinY = 0;
    int MaxX = SourceHeight;
    int MaxY = SourceWidth;
    
    tfit_bitmap SourceBitmap = TFIT_MakeBitmap(SourceAddress, SourceWidth, SourceHeight);
    
    void* Result = TFIT_Allocate(sizeof(tfit_u32)*SourceWidth*SourceHeight);
    if(Result)
    {
        int RotatedPitch = (int)(sizeof(tfit_u32)*SourceHeight);
        tfit_u8* Row = (tfit_u8*)Result;
        if(Direction == TFITRotationDir_Clockwise)
        {
            for(int Y = MinY;
                Y < MaxY;
                Y++)
            {
                tfit_u32* TexelPtr = (tfit_u32*)Row;
                for(int X = MinX;
                    X < MaxX;
                    X++)
                {
                    tfit_u32 SourceTexel =
                        *(tfit_u32*)TFIT_GET_BITMAP_PTR(SourceBitmap, Y, (MaxX - X - 1));
                    *TexelPtr++ = SourceTexel;
                }
                
                Row += RotatedPitch;
            }
        }
        else if(Direction == TFITRotationDir_Counterclockwise)
        {
            for(int Y = MinY;
                Y < MaxY;
                Y++)
            {
                tfit_u32* TexelPtr = (tfit_u32*)Row;
                for(int X = MinX;
                    X < MaxX;
                    X++)
                {
                    tfit_u32 SourceTexel =
                        *(tfit_u32*)TFIT_GET_BITMAP_PTR(SourceBitmap, (MaxY - Y - 1), X);
                    *TexelPtr++ = SourceTexel;
                }
                
                Row += RotatedPitch;
            }
        }
    }
    
    return(Result);
}

static void
TFITFree(void* Address)
{
    TFIT_Free(Address);
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
 furnished to do so, subject to the following conditions:
 
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
