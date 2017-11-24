// Copyright (C) 2017, Dominik Wodniok
// This software may be modified and distributed under the terms
// of the BSD 3-Clause license.
// See the LICENSE.txt file for details.

#ifndef EXAMPLE_H_INCLUDED
#define EXAMPLE_H_INCLUDED

/// \file   example.h
/// \author Dominik Wodniok
/// \date   2009

// std includes
#include <string>

// stl includes
#include <vector>

// dual mc builder vertex and quad definitions
#include "dualmc.h"

/// Example application for demonstrating the dual marching cubes builder.
class DualMCExample {
public:
    /// run example
    void run(int const argc, char** argv); 
    
private:

    /// Structure for the program options.
    struct AppOptions {
        std::string inputFile;
        int32_t dimX;
        int32_t dimY;
        int32_t dimZ;
        float isoValue;
        bool generateCaffeine;
        bool generateQuadSoup;
        std::string outputFile;
    };

    /// Parse program arguments.
    bool parseArgs(int const argc, char** argv, AppOptions & options);

    /// Generate an example volume for the dual mc builder.
    void generateCaffeine();
    
    /// Load volume from raw file.
    bool loadRawFile(std::string const & fileName, int32_t dimX, int32_t dimY, int32_t dimZ);

    /// Compute the iso surface for the specified iso value. Optionally generate
    /// a quad soup.
    void computeSurface(float const iso, bool const generateSoup = false);
    
    /// Write a Wavefront OBJ model for the extracted ISO surface.
    void writeOBJ(std::string const & fileName) const;
    
    /// Print program arguments.
    void printArgs() const;
    
    /// Print program help hint.
    void printHelpHint() const;
   
private:
    /// struct for volume data information
    struct Volume {
        // volume grid extents
        int32_t dimX;
        int32_t dimY;
        int32_t dimZ;
        // bit depth, should be 8 or 16
        int32_t bitDepth;
        /// volume data
        std::vector<uint8_t> data;
    };
       
    /// example volume
    Volume volume;
    
    /// Class for a volumetric sphere with gaussian fall-off.
    class RadialGaussian {
    public:
        /// Initialize with center coordinates and half density radius.
        RadialGaussian(float cX, float cY, float cZ, float variance);
        // evaluate the sphere function
        float eval(float x, float y, float z) const;
    private:
        // Coordinates of the sphere center.
        float cX;
        float cY;
        float cZ;
        // precomputed factors
        float normalization;
        float falloff;
        
    };

    /// array of vertices for the extracted surface
    std::vector<dualmc::Vertex> vertices;
    
    /// array of quad indices for the extracted surface
    std::vector<dualmc::Quad> quads;
};    

#endif // EXAMPLE_H_INCLUDED