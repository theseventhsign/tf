// ============================================================================
//  
//    tf_color_space v0.1a - Color interpolator by @twelvefifteen on GitHub
// 
// USAGE
// 
//    The user-facing functions can be found in the bottommost code section of
//    this file.
//    
//    float OutR, OutG, OutB;
//    TFCSInterpolateFloat(1.0f, 0.0f, 0.0f,
//                         0.5f,
//                         0.0f, 0.0f, 1.0f,
//                         TFCSColorSpace_HSV,
//                         &OutR, &OutG, &OutB);
//    // This call will interpolate halfway between red and blue in HSV.
//    // Colors for this call are represented using floating-point numbers on 
//    // the range [0.0f, 1.0f].
//    // The result can be found in the out parameters OutR, OutG, OutB.
// 
//    unsigned int Triplet;
//    TFCSInterpolateTriplet(0xFF0000, 0.5f, 0x0000FF,
//                           TFCSColorSpace_HSV,
//                           &Triplet);
//    // This call also interpolates halfway between red and blue in HSV.
//    // Colors for the call are represented with RGB triplets.
//    // The result can be found in the out parameter Triplet.
// 
// COLOR SPACES
// 
//    The set of supported color spaces are defined in tfcs_color_space below.
// 
//    Some color spaces (HSL, HSV, and CIELCh) use hue, which represents the
//    azimuth along what is effectively a color wheel. As a result, there are
//    two paths along which hue can be interpolated: a short one and a long
//    one. This library interpolates along the short path by default, but the
//    long path can be used instead by choosing the "reverse" flavors of these
//    color spaces.
// 
// LICENSE
// 
//    The license for this library can be found at the bottom of this file.
// 
// ============================================================================

#ifndef TF_COLOR_SPACE_H
#define TF_COLOR_SPACE_H

#include <math.h> // NOTE: atan2f, cbrtf, floorf, powf, sqrtf

#define TFCS_PI32 3.14159265359f

#ifdef _MSC_VER
typedef unsigned int tfcs_u32;
#else
#include <stdint.h>
typedef uint32_t tfcs_u32;
#endif
typedef float tfcs_f32;

#define TFCS_ARRAY_COUNT(Array) (sizeof((Array)) / sizeof((Array)[0]))

#define TFCS_MAXIMUM(A, B) ((A) > (B) ? (A) : (B))
#define TFCS_MINIMUM(A, B) ((A) < (B) ? (A) : (B))

#define TFCS_TO_DEGREES(Radians) ((Radians)*(180.0f / TFCS_PI32))
#define TFCS_TO_RADIANS(Degrees) ((Degrees)*(TFCS_PI32 / 180.0f))

typedef enum tfcs_color_space
{
    TFCSColorSpace_sRGB,
    TFCSColorSpace_LinearRGB,
    
    TFCSColorSpace_HSL,
    TFCSColorSpace_HSV,
    
    TFCSColorSpace_CIELAB,
    TFCSColorSpace_CIELCh,
    
    TFCSColorSpace_ReverseHSL,
    TFCSColorSpace_ReverseHSV,
    TFCSColorSpace_ReverseCIELCh,
    
	TFCSColorSpace_Count,
} tfcs_color_space;

typedef union tfcs_v3
{
    tfcs_f32 E[3];
    struct
    {
        tfcs_f32 x, y, z;
    } c; // NOTE: "Channels", "components", etc.
} tfcs_v3;

typedef struct tfcs_m3x3
{
	tfcs_v3 XAxis;
	tfcs_v3 YAxis;
	tfcs_v3 ZAxis;
} tfcs_m3x3;

static tfcs_v3
TFCS_V3(tfcs_f32 X, tfcs_f32 Y, tfcs_f32 Z)
{
	tfcs_v3 Result;
    Result.c.x = X;
    Result.c.y = Y;
    Result.c.z = Z;
    
	return(Result);
}

//
// NOTE: Scalar operations
//

static tfcs_f32
TFCS_ATan2(tfcs_f32 Y, tfcs_f32 X)
{
	tfcs_f32 Result = atan2f(Y, X);
    
	return(Result);
}

static tfcs_f32
TFCS_Clamp(tfcs_f32 Min, tfcs_f32 A, tfcs_f32 Max)
{
	tfcs_f32 Result = A;
	if (A < Min)
	{
		Result = Min;
	}
	else if (A > Max)
	{
		Result = Max;
	}
	return(Result);
}

static tfcs_f32
TFCS_Clamp01(tfcs_f32 A)
{
	tfcs_f32 Result = TFCS_Clamp(0.0f, A, 1.0f);
    
	return(Result);
}

static tfcs_f32
TFCS_Cos(tfcs_f32 A)
{
	tfcs_f32 Result = cosf(A);
    
	return(Result);
}

static tfcs_f32
TFCS_Cube(tfcs_f32 A)
{
	tfcs_f32 Result = (A*A*A);
    
	return(Result);
}

static tfcs_f32
TFCS_CubeRoot(tfcs_f32 A)
{
	tfcs_f32 Result = cbrtf(A);
    
	return(Result);
}

static tfcs_f32
TFCS_Floor(tfcs_f32 A)
{
	tfcs_f32 Result = floorf(A);
    
	return(Result);
}

static tfcs_f32
TFCS_Lerp(tfcs_f32 A, tfcs_f32 t, tfcs_f32 B)
{
	tfcs_f32 Result = (((1.0f - t)*A) + (t*B));
    
	return(Result);
}

static tfcs_f32
TFCS_MaxElement(tfcs_v3 A)
{
	tfcs_f32 Result = TFCS_MAXIMUM(A.c.x, TFCS_MAXIMUM(A.c.y, A.c.z));
    
	return(Result);
}

static tfcs_f32
TFCS_MinElement(tfcs_v3 A)
{
	tfcs_f32 Result = TFCS_MINIMUM(A.c.x, TFCS_MINIMUM(A.c.y, A.c.z));
    
	return(Result);
}

static tfcs_f32
TFCS_Pow(tfcs_f32 X, tfcs_f32 Y)
{
	tfcs_f32 Result = powf(X, Y);
    
	return(Result);
}

static tfcs_u32
TFCS_RoundToU32(tfcs_f32 A)
{
    tfcs_u32 Result = (tfcs_u32)roundf(A);
    
    return(Result); 
}

static tfcs_f32
TFCS_SafeRatioN(tfcs_f32 Dividend, tfcs_f32 Divisor, tfcs_f32 N)
{
	tfcs_f32 Result = N;
	if (Divisor != 0.0f)
	{
		Result = (Dividend / Divisor);
	}
	return(Result);
}

static tfcs_f32
TFCS_SafeRatio0(tfcs_f32 Dividend, tfcs_f32 Divisor)
{
	tfcs_f32 Result = TFCS_SafeRatioN(Dividend, Divisor, 0.0f);
    
	return(Result);
}

static tfcs_f32
TFCS_Sin(tfcs_f32 A)
{
	tfcs_f32 Result = sinf(A);
    
	return(Result);
}

static tfcs_f32
TFCS_Square(tfcs_f32 A)
{
	tfcs_f32 Result = (A*A);
    
	return(Result);
}

static tfcs_f32
TFCS_SquareRoot(tfcs_f32 A)
{
	tfcs_f32 Result = sqrtf(A);
    
	return(Result);
}

//
// NOTE: v3 operations
//

static tfcs_f32
TFCS_Dot(tfcs_v3 A, tfcs_v3 B)
{
	tfcs_f32 Result = ((A.c.x*B.c.x) +
                       (A.c.y*B.c.y) +
                       (A.c.z*B.c.z));
	return(Result);
}

static tfcs_v3
TFCS_UnpackTriplet(tfcs_u32 Triplet)
{
	tfcs_f32 Inv255 = (1.0f / 255.0f);
	tfcs_v3 Result;
	Result.c.x = (((Triplet >> 16) & 0xFF)*Inv255);
	Result.c.y = (((Triplet >> 8) & 0xFF)*Inv255);
	Result.c.z = (((Triplet >> 0) & 0xFF)*Inv255);
    
	return(Result);
}

static tfcs_u32
TFCS_PackV3(tfcs_v3 A)
{
    tfcs_u32 Result = ((TFCS_RoundToU32(A.c.x*255.0f) << 16) |
                       (TFCS_RoundToU32(A.c.y*255.0f) << 8) |
                       (TFCS_RoundToU32(A.c.z*255.0f) << 0));
    return(Result);
}

static tfcs_v3
TFCS_LerpV3(tfcs_v3 A, tfcs_f32 t, tfcs_v3 B)
{
	tfcs_v3 Result;
	Result.c.x = TFCS_Lerp(A.c.x, t, B.c.x);
	Result.c.y = TFCS_Lerp(A.c.y, t, B.c.y);
	Result.c.z = TFCS_Lerp(A.c.z, t, B.c.z);
    
	return(Result);
}

static tfcs_f32 TFCS_sRGBToLinearScalar(tfcs_f32 Channel);
static tfcs_v3
TFCS_sRGBToLinear(tfcs_v3 A)
{
	tfcs_v3 Result;
	Result.c.x = TFCS_sRGBToLinearScalar(A.c.x);
	Result.c.y = TFCS_sRGBToLinearScalar(A.c.y);
	Result.c.z = TFCS_sRGBToLinearScalar(A.c.z);
    
	return(Result);
}

static tfcs_f32 TFCS_LinearTosRGBScalar(tfcs_f32 Channel);
static tfcs_v3
TFCS_LinearTosRGB(tfcs_v3 A)
{
	tfcs_v3 Result;
	Result.c.x = TFCS_LinearTosRGBScalar(A.c.x);
	Result.c.y = TFCS_LinearTosRGBScalar(A.c.y);
	Result.c.z = TFCS_LinearTosRGBScalar(A.c.z);
    
	return(Result);
}

static tfcs_f32
TFCS_ShortestArcDistance(tfcs_v3 A, tfcs_v3 B, int EIndex)
{
	tfcs_f32 ChannelA = A.E[EIndex];
	tfcs_f32 ChannelB = B.E[EIndex];
    
	tfcs_f32 Result = (ChannelB - ChannelA);
    
	if ((ChannelB > ChannelA) &&
		((ChannelB - ChannelA) > 180.0f))
	{
		Result = (ChannelB - (ChannelA + 360.0f));
	}
	else if ((ChannelB < ChannelA) &&
             ((ChannelA - ChannelB) > 180.0f))
	{
		Result = (ChannelB + 360.0f - ChannelA);
	}
    
	return(Result);
}

static tfcs_f32
TFCS_LongestArcDistance(tfcs_v3 A, tfcs_v3 B, int EIndex)
{
	tfcs_f32 Result = 0.0f;
    
	tfcs_f32 ChannelA = A.E[EIndex];
	tfcs_f32 ChannelB = B.E[EIndex];
    
	if (ChannelB > ChannelA)
	{
		tfcs_f32 d = (ChannelB - ChannelA);
		if (d > 180.0f)
		{
			Result = d;
		}
		else
		{
			Result = -(360.0f - d);
		}
	}
	else
	{
		tfcs_f32 d = (ChannelA - ChannelB);
		if (d > 180.0f)
		{
			Result = -d;
		}
		else
		{
			Result = (360.0f - d);
		}
	}
    
	return(Result);
}

//
// NOTE: m3x3 operations
//

static tfcs_v3
TFCS_MultiplyM3X3(tfcs_m3x3 A, tfcs_v3 V)
{
	tfcs_v3 Result;
	Result.c.x = TFCS_Dot(TFCS_V3(A.XAxis.c.x, A.YAxis.c.x, A.ZAxis.c.x), V);
	Result.c.y = TFCS_Dot(TFCS_V3(A.XAxis.c.y, A.YAxis.c.y, A.ZAxis.c.y), V);
	Result.c.z = TFCS_Dot(TFCS_V3(A.XAxis.c.z, A.YAxis.c.z, A.ZAxis.c.z), V);
    
	return(Result);
}

//
// NOTE: Color space operations
// 

// NOTE: White point coords for Illuminant D65
#define TFCS_Xn 0.950470f
#define TFCS_Yn 1.0f
#define TFCS_Zn 1.088830f

static tfcs_f32
TFCS_sRGBToLinearScalar(tfcs_f32 C)
{
    C = TFCS_Clamp01(C);
    
    tfcs_f32 Result;
    if(C <= 0.04045f)
    {
        Result = (C / 12.92f);
    }
    else
    {
        Result = TFCS_Pow(((C + 0.055f) / 1.055f), 2.4f);
    }
    
    return(Result);
}

static tfcs_f32
TFCS_LinearTosRGBScalar(tfcs_f32 C)
{
    C = TFCS_Clamp01(C);
    
    tfcs_f32 Result;
    if(C <= 0.0031308f)
    {
        Result = (12.92f*C);
    }
    else
    {
        Result = ((1.055f*TFCS_Pow(C, (1.0f / 2.4f))) - 0.055f);
    }
    
    return(Result);
}

static tfcs_v3
TFCS_sRGBToHSL(tfcs_v3 A)
{
    tfcs_v3 Result = { 0 };
    
    tfcs_f32 MaxChannel = TFCS_MaxElement(A);
    tfcs_f32 MinChannel = TFCS_MinElement(A);
    tfcs_f32 Range = (MaxChannel - MinChannel);
    
    // NOTE: Hue
    
    if(MaxChannel == MinChannel)
    {
        Result.c.x = 0.0f;
    }
    else if(MaxChannel == A.c.x)
    {
        Result.c.x = (60.0f*(0.0f + ((A.c.y - A.c.z) / Range)));
    }
    else if(MaxChannel == A.c.y)
    {
        Result.c.x = (60.0f*(2.0f + ((A.c.z - A.c.x) / Range)));
    }
    else if(MaxChannel == A.c.z)
    {
        Result.c.x = (60.0f*(4.0f + ((A.c.x - A.c.y) / Range)));
    }
    if(Result.c.x < 0.0f)
    {
        Result.c.x += 360.0f;
    }
    
    // NOTE: Lightness
    
    tfcs_f32 Lightness = ((MaxChannel + MinChannel)*0.5f); 
    
    Result.c.z = Lightness;
    
    // NOTE: Saturation
    
    if(MaxChannel == 0.0f)
    {
        Result.c.y = 0.0f;
    }
    else if(MinChannel == 1.0f)
    {
        Result.c.y = 0.0f;
    }
    else
    {
        Result.c.y =
            TFCS_SafeRatio0((MaxChannel - Lightness), TFCS_MINIMUM(Lightness, (1.0f - Lightness)));
    }
    
    return(Result);
}

static tfcs_v3
TFCS_HSLTosRGB(tfcs_v3 A)
{
    tfcs_v3 Result = { 0 };
    
    if(A.c.y == 0.0f)
    {
        Result.c.x = A.c.z;
        Result.c.y = A.c.z;
        Result.c.z = A.c.z;
    }
    else
    {
        tfcs_f32 q = 0.0f;
        if(A.c.z < 0.5f)
        {
            q = (A.c.z*(1.0f + A.c.y));
        }
        else if (A.c.z >= 0.5f)
        {
            q = (A.c.z + A.c.y - (A.c.z*A.c.y));
        }
        tfcs_f32 p = (2.0f*A.c.z - q);
        tfcs_f32 hSubk = (A.c.x / 360.0f);
        
        tfcs_f32 tSubC[3];
        tSubC[0] = (hSubk + (1.0f / 3.0f)); 
        tSubC[1] = hSubk;
        tSubC[2] = (hSubk - (1.0f / 3.0f)); 
        
        for(int tIndex = 0;
            tIndex < TFCS_ARRAY_COUNT(tSubC);
            tIndex++)
        {
            tfcs_f32* t = (tSubC + tIndex);
            if(*t < 0.0f)
            {
                *t += 1.0f;
            }
            if(*t > 1.0f)
            {
                *t -= 1.0f;
            }
        }
        
        for(int EIndex = 0;
            EIndex < TFCS_ARRAY_COUNT(Result.E);
            EIndex++)
        {
            tfcs_f32* Channel = (Result.E + EIndex);
            tfcs_f32* t = (tSubC + EIndex);
            
            if(*t < (1.0f / 6.0f))
            {
                *Channel = (p + ((q - p)*6.0f*(*t)));
            }
            else if(((1.0f / 6.0f) <= *t) && (*t < 0.5f))
            {
                *Channel = q;
            }
            else if((0.5f <= *t) && (*t < (2.0 / 3.0f)))
            {
                *Channel = (p + ((q - p)*6.0f*((2.0f / 3.0f) - *t)));
            }
            else
            {
                *Channel = p;
            }
        }
    }
    
    return(Result);
}

static tfcs_v3
TFCS_sRGBToHSV(tfcs_v3 A)
{
    tfcs_v3 Result = { 0 };
    
    tfcs_f32 MaxChannel = TFCS_MaxElement(A);
    tfcs_f32 MinChannel = TFCS_MinElement(A);
    tfcs_f32 Range = (MaxChannel - MinChannel);
    
    // NOTE: Hue
    
    if(MaxChannel == MinChannel)
    {
        Result.c.x = 0.0f;
    }
    else if(MaxChannel == A.c.x)
    {
        Result.c.x = (60.0f*(0.0f + ((A.c.y - A.c.z) / Range)));
    }
    else if(MaxChannel == A.c.y)
    {
        Result.c.x = (60.0f*(2.0f + ((A.c.z - A.c.x) / Range)));
    }
    else if(MaxChannel == A.c.z)
    {
        Result.c.x = (60.0f*(4.0f + ((A.c.x - A.c.y) / Range)));
    }
    if(Result.c.x < 0.0f)
    {
        Result.c.x += 360.0f;
    }
    
    // NOTE: Saturation
    
    if(MaxChannel == 0.0f)
    {
        Result.c.y = 0.0f;
    }
    else
    {
        Result.c.y = (Range / MaxChannel);
    }
    
    // NOTE: Value
    
    Result.c.z = MaxChannel;
    
    return(Result);
}

static tfcs_v3
TFCS_HSVTosRGB(tfcs_v3 A)
{
    tfcs_v3 Result = { 0 };
    
    if(A.c.y == 0.0f)
    {
        Result.c.x = A.c.z;
        Result.c.y = A.c.z;
        Result.c.z = A.c.z;
    }
    else
    {
        if(A.c.x == 360.0f)
        {
            A.c.x = 0.0f;
        }
        if(A.c.x > 360.0f)
        {
            A.c.x -= 360.0f;
        }
        if(A.c.x < 0.0f)
        {
            A.c.x += 360.0f;
        }
        
        A.c.x /= 60.0f;
        int I = (int)TFCS_Floor(A.c.x);
        tfcs_f32 F = (A.c.x - I);
        tfcs_f32 P = (A.c.z*(1.0f - A.c.y));
        tfcs_f32 Q = (A.c.z*(1.0f - (A.c.y*F)));
        tfcs_f32 T = (A.c.z*(1.0f - (A.c.y*(1.0f - F))));
        switch(I)
        {
            case 0:
            {
                Result.c.x = A.c.z;
                Result.c.y = T;
                Result.c.z = P;
            } break;
            
            case 1:
            {
                Result.c.x = Q;
                Result.c.y = A.c.z;
                Result.c.z = P;
            } break;
            
            case 2:
            {
                Result.c.x = P;
                Result.c.y = A.c.z;
                Result.c.z = T;
            } break;
            
            case 3:
            {
                Result.c.x = P;
                Result.c.y = Q;
                Result.c.z = A.c.z;
            } break;
            
            case 4:
            {
                Result.c.x = T;
                Result.c.y = P;
                Result.c.z = A.c.z;
            } break;
            
            case 5:
            {
                Result.c.x = A.c.z;
                Result.c.y = P;
                Result.c.z = Q;
            } break;
        }
    }
    
    return(Result);
}

static tfcs_v3
LinearToCIEXYZ(tfcs_v3 A)
{
    tfcs_m3x3 Matrix;
    Matrix.XAxis = TFCS_V3(0.4124564f, 0.2126729f, 0.0193339f);
    Matrix.YAxis = TFCS_V3(0.3575761f, 0.7151522f, 0.1191920f);
    Matrix.ZAxis = TFCS_V3(0.1804375f, 0.0721750f, 0.9503041f);
    
    tfcs_v3 Result = TFCS_MultiplyM3X3(Matrix, A);
    
    return(Result);
}

static tfcs_f32
TFCS_Ft(tfcs_f32 t)
{
    tfcs_f32 Result;
    
    tfcs_f32 Sigma = (6.0f / 29.0f);
    if(t > TFCS_Cube(Sigma))
    {
        Result = TFCS_CubeRoot(t);
    }
    else
    {
        Result = ((t / (3*TFCS_Square(Sigma))) + (4.0f / 29.0f));
    }
    
    return(Result);
}

static tfcs_v3
CIEXYZToCIELAB(tfcs_v3 A)
{
    tfcs_v3 Result;
    Result.c.x = (116.0f*TFCS_Ft(A.c.y / TFCS_Yn) - 16.0f);
    Result.c.y = (500.0f*(TFCS_Ft(A.c.x / TFCS_Xn) - TFCS_Ft(A.c.y / TFCS_Yn)));
    Result.c.z = (200.0f*(TFCS_Ft(A.c.y / TFCS_Yn) - TFCS_Ft(A.c.z / TFCS_Zn)));
    
    return(Result);
}

static tfcs_f32
TFCS_InvFt(tfcs_f32 t)
{
    tfcs_f32 Result;
    
    tfcs_f32 Sigma = (6.0f / 29.0f);
    if(t > Sigma)
    {
        Result = TFCS_Cube(t);
    }
    else
    {
        Result = ((3.0f*TFCS_Square(Sigma))*(t - (4.0f / 29.0f)));
    }
    
    return(Result);
}

static tfcs_v3
TFCS_CIELABToCIEXYZ(tfcs_v3 A)
{
    tfcs_v3 Result;
    Result.c.x = (TFCS_Xn*TFCS_InvFt(((A.c.x + 16.0f) / 116.0f) + (A.c.y / 500.0f)));
    Result.c.y = (TFCS_Yn*TFCS_InvFt((A.c.x + 16.0f) / 116.0f));
    Result.c.z = (TFCS_Zn*TFCS_InvFt(((A.c.x + 16.0f) / 116.0f) - (A.c.z / 200.0f)));
    
    return(Result);
}

static tfcs_v3
CIELABToCIELCh(tfcs_v3 A)
{
    tfcs_v3 Result;
    Result.c.x = A.c.x;
    Result.c.y = TFCS_SquareRoot(TFCS_Square(A.c.y) + TFCS_Square(A.c.z));
    tfcs_f32 C = TFCS_TO_DEGREES(TFCS_ATan2(A.c.z, A.c.y));
    if(C >= 0.0f)
    {
        Result.c.z = C;
    }
    else
    {
        Result.c.z = (C + 360.0f);
    }
    
    return(Result);
}

static tfcs_v3
TFCS_CIELChToCIELAB(tfcs_v3 A)
{
    tfcs_v3 Result;
    Result.c.x = A.c.x;
    Result.c.y = (A.c.y*TFCS_Cos(TFCS_TO_RADIANS(A.c.z)));
    Result.c.z = (A.c.y*TFCS_Sin(TFCS_TO_RADIANS(A.c.z)));
    
    return(Result);
}

static tfcs_v3
TFCS_CIEXYZToLinear(tfcs_v3 A)
{
    tfcs_m3x3 Matrix;
    Matrix.XAxis = TFCS_V3(3.2404542f, -0.9692660f, 0.0556434f);
    Matrix.YAxis = TFCS_V3(-1.5371385f, 1.8760108f, -0.2040259f);
    Matrix.ZAxis = TFCS_V3(-0.4985314f, 0.0415560f, 1.0572252f);
    
    tfcs_v3 Result = TFCS_MultiplyM3X3(Matrix, A);
    
    return(Result);
}

static tfcs_v3
TFCS_InterpolateV3(tfcs_v3 A, tfcs_f32 t, tfcs_v3 B, tfcs_color_space ColorSpace)
{
    tfcs_v3 Result = { 0 };
    
    switch(ColorSpace)
    {
        case TFCSColorSpace_sRGB:
        {
            Result = TFCS_LerpV3(A, t, B);
        } break;
        
        case TFCSColorSpace_LinearRGB:
        {
            A = TFCS_sRGBToLinear(A);
            B = TFCS_sRGBToLinear(B);
            
            Result = TFCS_LerpV3(A, t, B);
            Result = TFCS_LinearTosRGB(Result);
        } break;
        
        case TFCSColorSpace_HSL:
        {
            A = TFCS_sRGBToHSL(A);
            B = TFCS_sRGBToHSL(B);
            
            tfcs_v3 ResultHSL;
            tfcs_f32 dH = TFCS_ShortestArcDistance(A, B, 0);
            tfcs_f32 Hue = (A.c.x + (t*dH));
            ResultHSL.c.x = Hue;
            ResultHSL.c.y = TFCS_Lerp(A.c.y, t, B.c.y);
            ResultHSL.c.z = TFCS_Lerp(A.c.z, t, B.c.z);
            
            Result = TFCS_HSLTosRGB(ResultHSL);
        } break;
        
        case TFCSColorSpace_ReverseHSL:
        {
            A = TFCS_sRGBToHSL(A);
            B = TFCS_sRGBToHSL(B);
            
            tfcs_v3 ResultHSL;
            tfcs_f32 dH = TFCS_LongestArcDistance(A, B, 0);
            tfcs_f32 Hue = (A.c.x + (t*dH));
            ResultHSL.c.x = Hue;
            ResultHSL.c.y = TFCS_Lerp(A.c.y, t, B.c.y);
            ResultHSL.c.z = TFCS_Lerp(A.c.z, t, B.c.z);
            
            Result = TFCS_HSLTosRGB(ResultHSL);
        } break;
        
        case TFCSColorSpace_HSV:
        {
            A = TFCS_sRGBToHSV(A);
            B = TFCS_sRGBToHSV(B);
            
            tfcs_v3 ResultHSV;
            tfcs_f32 dH = TFCS_ShortestArcDistance(A, B, 0);
            tfcs_f32 Hue = (A.c.x + (t*dH));
            ResultHSV.c.x = Hue;
            ResultHSV.c.y = TFCS_Lerp(A.c.y, t, B.c.y);
            ResultHSV.c.z = TFCS_Lerp(A.c.z, t, B.c.z);
            
            Result = TFCS_HSVTosRGB(ResultHSV);
        } break;
        
        case TFCSColorSpace_ReverseHSV:
        {
            A = TFCS_sRGBToHSV(A);
            B = TFCS_sRGBToHSV(B);
            
            tfcs_v3 ResultHSV;
            tfcs_f32 dH = TFCS_LongestArcDistance(A, B, 0);
            tfcs_f32 Hue = (A.c.x + (t*dH));
            ResultHSV.c.x = Hue;
            ResultHSV.c.y = TFCS_Lerp(A.c.y, t, B.c.y);
            ResultHSV.c.z = TFCS_Lerp(A.c.z, t, B.c.z);
            
            Result = TFCS_HSVTosRGB(ResultHSV);
        } break;
        
        case TFCSColorSpace_CIELAB:
        {
            A = TFCS_sRGBToLinear(A);
            B = TFCS_sRGBToLinear(B);
            A = LinearToCIEXYZ(A);
            B = LinearToCIEXYZ(B);
            A = CIEXYZToCIELAB(A);
            B = CIEXYZToCIELAB(B);
            
            tfcs_v3 ResultCIELAB = TFCS_LerpV3(A, t, B);
            tfcs_v3 ResultCIEXYZ = TFCS_CIELABToCIEXYZ(ResultCIELAB);
            tfcs_v3 ResultLinearRGB = TFCS_CIEXYZToLinear(ResultCIEXYZ);
            Result = TFCS_LinearTosRGB(ResultLinearRGB);
        } break;
        
        case TFCSColorSpace_CIELCh:
        {
            A = TFCS_sRGBToLinear(A);
            B = TFCS_sRGBToLinear(B);
            A = LinearToCIEXYZ(A);
            B = LinearToCIEXYZ(B);
            A = CIEXYZToCIELAB(A);
            B = CIEXYZToCIELAB(B);
            A = CIELABToCIELCh(A);
            B = CIELABToCIELCh(B);
            
            tfcs_v3 ResultCIELCh;
            ResultCIELCh.c.x = TFCS_Lerp(A.c.x, t, B.c.x);
            ResultCIELCh.c.y = TFCS_Lerp(A.c.y, t, B.c.y);
            tfcs_f32 dH = TFCS_ShortestArcDistance(A, B, 2);
            ResultCIELCh.c.z = (A.c.z + (t*dH));
            
            tfcs_v3 ResultCIELAB = TFCS_CIELChToCIELAB(ResultCIELCh);
            tfcs_v3 ResultCIEXYZ = TFCS_CIELABToCIEXYZ(ResultCIELAB);
            tfcs_v3 ResultLinearRGB = TFCS_CIEXYZToLinear(ResultCIEXYZ);
            Result = TFCS_LinearTosRGB(ResultLinearRGB);
        } break;
        
        case TFCSColorSpace_ReverseCIELCh:
        {
            A = TFCS_sRGBToLinear(A);
            B = TFCS_sRGBToLinear(B);
            A = LinearToCIEXYZ(A);
            B = LinearToCIEXYZ(B);
            A = CIEXYZToCIELAB(A);
            B = CIEXYZToCIELAB(B);
            A = CIELABToCIELCh(A);
            B = CIELABToCIELCh(B);
            
            tfcs_v3 ResultCIELCh;
            ResultCIELCh.c.x = TFCS_Lerp(A.c.x, t, B.c.x);
            ResultCIELCh.c.y = TFCS_Lerp(A.c.y, t, B.c.y);
            tfcs_f32 dH = TFCS_LongestArcDistance(A, B, 2);
            ResultCIELCh.c.z = (A.c.z + (t*dH));
            
            tfcs_v3 ResultCIELAB = TFCS_CIELChToCIELAB(ResultCIELCh);
            tfcs_v3 ResultCIEXYZ = TFCS_CIELABToCIEXYZ(ResultCIELAB);
            tfcs_v3 ResultLinearRGB = TFCS_CIEXYZToLinear(ResultCIEXYZ);
            Result = TFCS_LinearTosRGB(ResultLinearRGB);
        } break;
        
        case TFCSColorSpace_Count:
        {
            
        } break;
    }
    
    return(Result);
}

//
// NOTE: User-facing API
//

static void
TFCSInterpolateFloat(tfcs_f32 AR, tfcs_f32 AG, tfcs_f32 AB,
                     tfcs_f32 t,
                     tfcs_f32 BR, tfcs_f32 BG, tfcs_f32 BB,
                     tfcs_color_space ColorSpace,
                     tfcs_f32* OutR, tfcs_f32* OutG, tfcs_f32* OutB)
{
    t = TFCS_Clamp01(t);
    
    tfcs_v3 A = TFCS_V3(AR, AG, AB);
    tfcs_v3 B = TFCS_V3(BR, BG, BB);
    
    tfcs_v3 ResultV3 = TFCS_InterpolateV3(A, t, B, ColorSpace);
    
    *OutR = ResultV3.c.x;
    *OutG = ResultV3.c.y;
    *OutB = ResultV3.c.z;
}

static void
TFCSInterpolateTriplet(tfcs_u32 TripletA, tfcs_f32 t, tfcs_u32 TripletB,
                       tfcs_color_space ColorSpace,
                       tfcs_u32* OutTriplet)
{
    t = TFCS_Clamp01(t);
    
    tfcs_v3 A = TFCS_UnpackTriplet(TripletA);
    tfcs_v3 B = TFCS_UnpackTriplet(TripletB);
    
    tfcs_v3 ResultV3 = TFCS_InterpolateV3(A, t, B, ColorSpace);
    
    *OutTriplet = TFCS_PackV3(ResultV3);
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
