// Copyright (C) 2017, Dominik Wodniok
// This software may be modified and distributed under the terms
// of the BSD 3-Clause license.
// See the LICENSE.txt file for details.

/// \file   example.cpp
/// \author Dominik Wodniok
/// \date   2009

/*

Modified by Adrian Bowyer to add loading of simple 3D tensor files of floats and to output ASCII STL files
in addition to the original Wavefront OBJ files.

TODO: Output binary STL as well.

Dr Adrian Bowyer
RepRap Ltd
https://reprapltd.com
contact@reprapltd.com

12 August 2019

The tensor file format is a list of three integers then floats all separated by spaces:

xDimension yDimension zDimension minValue maxValue { (xDimension*yDimension*zDimension) function values } EoF

The first three integers are the size of the tensor.  The second two floats are the maximum and minimum values of all the
numbers in the tensor, then there is a list of the float values incrementing x fastest, then y, then z.

Example command line:

$ ./dmc -tensor data/testCylinder.tns -iso 0.5 -out data/testCylinder.stl

*/

// C libs
#include <cmath>
#include <cstdlib>
#include <cstring>

// std libs
#include <chrono>
#include <fstream>
#include <iostream>
#include <limits>

// stl
#include <vector>

#include<algorithm>

// dual mc builder
#include "dualmc.h"

// main include
#include "example.h"

using std::chrono::high_resolution_clock;
using std::chrono::duration;
using std::chrono::duration_cast;

bool debug = true;

//------------------------------------------------------------------------------

void DualMCExample::run(int const argc, char** argv) {
    // parse program options
    AppOptions options;
    if(!parseArgs(argc,argv,options)) {
        return;
    }
    
    // load raw file or generate example volume dataset
    if(options.generateCaffeine) {
        generateCaffeine();
    } else if(!options.inputFile.empty()) {
        if(options.readTensor)
	{
	  if(!loadTensor(options.inputFile)) {
            return;
          }
        } else if(!loadRawFile(options.inputFile, options.dimX, options.dimY, options.dimZ)) {
            return;
        }
    } else {
        std::cerr << "No input specified" << std::endl;
        printHelpHint();
        return;
    }
    
    // compute ISO surface
    if(debug)
        std::cout << "Computing iso-surface." << std::endl;
    computeSurface(options.isoValue,options.generateQuadSoup,options.generateManifold);

    std::string extn = options.outputFile.substr(options.outputFile.size()-4,options.outputFile.size()-1);
    std::transform(extn.begin(), extn.end(), extn.begin(), ::tolower);
    
    // write output file
    if(debug)
        std::cout << "Writing output file." << std::endl;

    if(extn.compare(".obj") == 0)
    	writeOBJ(options.outputFile);
    else if(extn.compare(".stl") == 0)
	writeSTL(options.outputFile);
    else
	std::cerr << "Output file is neither .obj nor .stl: " << options.outputFile << std::endl;
}

//------------------------------------------------------------------------------

bool DualMCExample::parseArgs(int const argc, char** argv, AppOptions & options) {
    // set default values
    options.inputFile.assign("");
    options.dimX = -1;
    options.dimY = -1;
    options.dimZ = -1;
    options.isoValue = 0.5f;
    options.generateCaffeine = false;
    options.generateQuadSoup = false;
    options.generateManifold = false;
    options.readTensor = false;
    options.outputFile.assign("surface.obj");
    
    // parse arguments
    for(int currentArg = 1; currentArg < argc; ++currentArg) {
        if(strcmp(argv[currentArg],"-soup") == 0) {
            options.generateQuadSoup = true;
        } else if(strcmp(argv[currentArg],"-caffeine") == 0) {
            options.generateCaffeine = true;
        } else if(strcmp(argv[currentArg],"-manifold") == 0) {
            options.generateManifold = true;
        } else if(strcmp(argv[currentArg],"-iso") == 0) {
            if(currentArg+1 == argc) {
                std::cerr << "Iso value missing" << std::endl;
                return false;
            }
            // Read the iso value and clamp it to [0,1].
            // Invalid values are set to 0.
            options.isoValue = atof(argv[currentArg+1]);
            if(options.isoValue > 1.0f)
                options.isoValue = 1.0f;
            else if(options.isoValue < 0.0f || options.isoValue != options.isoValue)
                options.isoValue = 0.0f;
            ++currentArg;
        } else if(strcmp(argv[currentArg],"-out") == 0) {
            if(currentArg+1 == argc) {
                std::cerr << "Output filename missing" << std::endl;
                return false;
            }
            options.outputFile.assign(argv[currentArg+1]);
            ++currentArg;
        } else if(strcmp(argv[currentArg],"-raw") == 0) {
            if(currentArg+4 >= argc) {
                std::cerr << "Not enough arguments for raw file" << std::endl;
                return false;
            }
            options.inputFile.assign(argv[currentArg+1]);
            options.dimX = atoi(argv[currentArg+2]);
            options.dimY = atoi(argv[currentArg+3]);
            options.dimZ = atoi(argv[currentArg+4]);
            currentArg += 4;
        }else if(strcmp(argv[currentArg],"-tensor") == 0) {
            if(currentArg+1 >= argc) {
                std::cerr << "Not enough arguments for tensor file" << std::endl;
                return false;
            }
            options.inputFile.assign(argv[currentArg+1]);
            currentArg += 1;
	    options.readTensor = true;
        } else if(strcmp(argv[currentArg],"-help") == 0) {
            printArgs();
            return false;
        } else {
            std::cerr << "Unknown argument: " << argv[currentArg] << std::endl;
            printHelpHint();
            return false;
        }
    }
    return true;
}

//------------------------------------------------------------------------------

void DualMCExample::printArgs() const {
    std::cout << "Usage: dmc ARGS" << std::endl;
    std::cout << " -help              print this help" << std::endl;
    std::cout << " -raw FILE X Y Z    specify raw file with dimensions" << std::endl;
    std::cout << " -tensor FILE       specify tensor file" << std::endl;
    std::cout << " -caffeine          generate built-in caffeine molecule" << std::endl;
    std::cout << " -manifold          use Manifold Dual Marching Cubes algorithm (Rephael Wenger)" << std::endl;
    std::cout << " -iso X             specify iso value X in [0,1]. DEFAULT: 0.5" << std::endl;
    std::cout << " -out FILE          specify output file name. Extension (.obj or .stl) determines file type. DEFAULT: surface.obj" << std::endl;
    std::cout << " -soup              generate a quad soup (no vertex sharing)" << std::endl;
}

//------------------------------------------------------------------------------

void DualMCExample::printHelpHint() const {
    std::cout << "Try: dmc -help" << std::endl;
}

//------------------------------------------------------------------------------

void DualMCExample::computeSurface(float const iso, bool const generateSoup, bool const generateManifold) {
    std::cout << "Computing surface" << std::endl;
    
    // measure extraction time
    high_resolution_clock::time_point const startTime = high_resolution_clock::now();
    
    // construct iso surface
    if(volume.bitDepth == 8) {
        dualmc::DualMC<uint8_t> builder;
        builder.build(&volume.data.front(), volume.dimX, volume.dimY, volume.dimZ,
            iso * std::numeric_limits<uint8_t>::max(), generateManifold, generateSoup, vertices, quads);
    } else if(volume.bitDepth == 16) {
        dualmc::DualMC<uint16_t> builder;
        builder.build((uint16_t const*)&volume.data.front(), volume.dimX, volume.dimY, volume.dimZ,
            iso * std::numeric_limits<uint16_t>::max(), generateManifold, generateSoup, vertices, quads);
    } else {
        std::cerr << "Invalid volume bit depth" << std::endl;
        return;
    }
        
    high_resolution_clock::time_point const endTime = high_resolution_clock::now();
    duration<double> const diffTime = duration_cast<duration<double>>(endTime - startTime);
    double const extractionTime = diffTime.count();
    
    std::cout << "Extraction time: " << 1000*extractionTime << "ms" << std::endl;
}

//------------------------------------------------------------------------------

void DualMCExample::generateCaffeine() {
    std::cout << "Generating caffeine volume" << std::endl;
    
    // initialize volume dimensions and memory
    volume.dimX = 128;
    volume.dimY = 128;
    volume.dimZ = 128;
    size_t const numDataPoints = volume.dimX * volume.dimY * volume.dimZ;
    volume.data.resize(numDataPoints*2);
    volume.bitDepth = 16;
    
    float invDimX = 1.0f / (volume.dimX-1);
    float invDimY = 1.0f / (volume.dimY-1);
    float invDimZ = 1.0f / (volume.dimZ-1);
    
    // create caffeine molecule
    // 3D structure from https://pubchem.ncbi.nlm.nih.gov/compound/caffeine#section=Top
    
    // caffeine scale
    float constexpr s = 1.0f/10.0f;
    // caffeine offset
    float constexpr oX = 0.5f;
    float constexpr oY = 0.5f;
    float constexpr oZ = 0.5f;
    // atom scale scale
    //float constexpr as = 0.001f/70.0f/70.0f;
    float constexpr as = 0.025*0.025/70.0f/70.0f;
    // atom scales
    float const atomScales[] = {25*25*as,70*70*as,65*65*as,60*60*as};
    enum ElementType {HYDROGEN=0,CARBON=1,NITROGEN=2,OXYGEN=3};
    
    // approximate electron density with radial Gaussians.
    std::vector<RadialGaussian> atoms;
    atoms.reserve(24);
    // 1 hydrogen, 6 carbon, 7 nitrogen, 8 oxygen
    atoms.emplace_back(   0.47 * s + oX,  2.5688 * s + oY,  0.0006 * s + oZ,atomScales[OXYGEN]); // 8
    atoms.emplace_back(-3.1271 * s + oX, -0.4436 * s + oY, -0.0003 * s + oZ,atomScales[OXYGEN]); // 8
    atoms.emplace_back(-0.9686 * s + oX, -1.3125 * s + oY,       0 * s + oZ,atomScales[NITROGEN]); // 7
    atoms.emplace_back( 2.2182 * s + oX,  0.1412 * s + oY, -0.0003 * s + oZ,atomScales[NITROGEN]); // 7
    atoms.emplace_back(-1.3477 * s + oX,  1.0797 * s + oY, -0.0001 * s + oZ,atomScales[NITROGEN]); // 7
    atoms.emplace_back( 1.4119 * s + oX, -1.9372 * s + oY,  0.0002 * s + oZ,atomScales[NITROGEN]); // 7
    atoms.emplace_back( 0.8579 * s + oX,  0.2592 * s + oY, -0.0008 * s + oZ,atomScales[CARBON]); // 6
    atoms.emplace_back( 0.3897 * s + oX, -1.0264 * s + oY, -0.0004 * s + oZ,atomScales[CARBON]); // 6
    atoms.emplace_back(-1.9061 * s + oX, -0.2495 * s + oY, -0.0004 * s + oZ,atomScales[CARBON]); // 6
    atoms.emplace_back( 0.0307 * s + oX,   1.422 * s + oY, -0.0006 * s + oZ,atomScales[CARBON]); // 6
    atoms.emplace_back( 2.5032 * s + oX, -1.1998 * s + oY,  0.0003 * s + oZ,atomScales[CARBON]); // 6
    atoms.emplace_back(-1.4276 * s + oX, -2.6960 * s + oY,  0.0008 * s + oZ,atomScales[CARBON]); // 6
    atoms.emplace_back( 3.1926 * s + oX,  1.2061 * s + oY,  0.0003 * s + oZ,atomScales[CARBON]); // 6
    atoms.emplace_back(-2.2969 * s + oX,  2.1881 * s + oY,  0.0007 * s + oZ,atomScales[CARBON]); // 6
    atoms.emplace_back( 3.5163 * s + oX, -1.5787 * s + oY,  0.0008 * s + oZ,atomScales[HYDROGEN]); // 1
    atoms.emplace_back(-1.0451 * s + oX, -3.1973 * s + oY, -0.8937 * s + oZ,atomScales[HYDROGEN]); // 1
    atoms.emplace_back(-2.5186 * s + oX, -2.7596 * s + oY,  0.0011 * s + oZ,atomScales[HYDROGEN]); // 1
    atoms.emplace_back(-1.0447 * s + oX, -3.1963 * s + oY,  0.8957 * s + oZ,atomScales[HYDROGEN]); // 1
    atoms.emplace_back( 4.1992 * s + oX,  0.7801 * s + oY,  0.0002 * s + oZ,atomScales[HYDROGEN]); // 1
    atoms.emplace_back( 3.0468 * s + oX,  1.8092 * s + oY, -0.8992 * s + oZ,atomScales[HYDROGEN]); // 1
    atoms.emplace_back( 3.0466 * s + oX,  1.8083 * s + oY,  0.9004 * s + oZ,atomScales[HYDROGEN]); // 1
    atoms.emplace_back(-1.8087 * s + oX,  3.1651 * s + oY, -0.0003 * s + oZ,atomScales[HYDROGEN]); // 1
    atoms.emplace_back(-2.9322 * s + oX,  2.1027 * s + oY,  0.8881 * s + oZ,atomScales[HYDROGEN]); // 1
    atoms.emplace_back(-2.9346 * s + oX,  2.1021 * s + oY, -0.8849 * s + oZ,atomScales[HYDROGEN]); // 1
    
    uint16_t * data16Bit = (uint16_t*)&volume.data.front();
    
    // scale for density field
    float constexpr postDensityScale = 2.5f;
    
    // volume write position
    int32_t p = 0;
    // iterate all voxels
    // compute canoncical [0,1]^3 volume coordinates for density evaluation
    for(int32_t z = 0; z < volume.dimZ; ++z) {
        float const nZ = float(z) * invDimZ;
        for(int32_t y = 0; y < volume.dimY; ++y) {
            float const nY = float(y) * invDimY;
            for(int32_t x = 0; x < volume.dimX; ++x, ++p) {
                float const nX = float(x) * invDimX;
                float rho = 0.0f;
                // compute sum of electron densities
                for(auto const & a : atoms) {
                    rho += a.eval(nX,nY,nZ);
                }
                rho *= postDensityScale;
                if(rho > 1.0f)
                    rho = 1.0f;
                data16Bit[p] = rho * std::numeric_limits<uint16_t>::max();
            }
        }
    }
}

// Load a 3D tensor of float values and convert them to 16 bit densities.
 
bool DualMCExample::loadTensor(std::string const & fileName) {
    std::cout << "Loading tensor file " << fileName << std::endl;

    std::ifstream file(fileName, std::ifstream::in);
    if(!file) {
        std::cerr << "Unable to open file '" << fileName << "'" << std::endl;
        return false;
    }
    
    // initialize volume dimensions and memory
    file >> volume.dimX;
    file >> volume.dimY;
    file >> volume.dimZ;
    size_t const numDataPoints = volume.dimX * volume.dimY * volume.dimZ;
    volume.data.resize(numDataPoints*2);
    volume.bitDepth = 16;

    float minValue, maxValue, scale;
    file >> minValue;
    file >> maxValue;
    scale = 1.0/(maxValue - minValue);

    uint16_t * data16Bit = (uint16_t*)&volume.data.front();
    
    // volume write position
    int32_t p = 0;
    float rho;

    for(int32_t z = 0; z < volume.dimZ; ++z) {
        //float const nZ = float(z) * invDimZ;
        for(int32_t y = 0; y < volume.dimY; ++y) {
            //float const nY = float(y) * invDimY;
            for(int32_t x = 0; x < volume.dimX; ++x, ++p) {
                file >> rho;
		if(rho < minValue || rho > maxValue)
			std::cerr << "Value numbered " << p << " in the tensor is outside the range: " << rho << std::endl;
		rho = scale*(rho - minValue);
                data16Bit[p] = rho * std::numeric_limits<uint16_t>::max();
            }
        }
    }

    if(debug)
    	std::cout << "Number of values in the tensor is: " << p << std::endl;

    if(!file) {
        std::cerr << "Error while reading file" << std::endl;
        return false;
    } else
	file.close();

    return true;
}


//------------------------------------------------------------------------------

bool DualMCExample::loadRawFile(std::string const & fileName, int32_t dimX, int32_t dimY, int32_t dimZ) {
    // check provided dimensions
    if(dimX < 1 || dimY < 1 || dimZ < 1) {
        std::cerr << "Invalid RAW file dimensions specified" << std::endl;
        return false;
    }
    
    // open raw file
    std::ifstream file(fileName, std::ifstream::binary);
    if(!file) {
        std::cerr << "Unable to open file '" << fileName << "'" << std::endl;
        return false;
    }
    
    // check consistency of file size and volume dimensions
    size_t const expectedFileSize = size_t(dimX) * size_t(dimY) * size_t(dimZ);
    file.seekg (0, file.end);
    size_t const fileSize = file.tellg();
    file.seekg (0, file.beg);
    
    if(expectedFileSize != fileSize) {
        if(expectedFileSize * 2 == fileSize) {
            std::cout << "Assuming 16-bit RAW file" << std::endl;
            volume.bitDepth = 16;
        } else {
            std::cerr << "File size inconsistent with specified dimensions" << std::endl;
            return false;
        }
    } else {
        volume.bitDepth = 8;
    }

    //
    if(expectedFileSize >= 0xffffffffu) {
        std::cerr << "Too many voxels. Please improve the dual mc implementation." << std::endl;
        return false;
    }
    
    // initialize volume dimensions and memory
    volume.dimX = dimX;
    volume.dimY = dimY;
    volume.dimZ = dimZ;
    volume.data.resize(fileSize);
    
    // read data
    file.read((char*)&volume.data[0], fileSize);
    
    if(!file) {
        std::cerr << "Error while reading file" << std::endl;
        return false;
    }
    
    return true;
}

//------------------------------------------------------------------------------

void DualMCExample::writeOBJ(std::string const & fileName) const {
    std::cout << "Writing OBJ file " << fileName << std::endl;
    // check if we actually have an ISO surface
    if(vertices.size () == 0 || quads.size() == 0) {
        std::cout << "No ISO surface generated. Skipping OBJ generation." << std::endl;
        return;
    }
    
    // open output file
    std::ofstream file(fileName);
    if(!file) {
        std::cout << "Error opening output file" << std::endl;
        return;
    }
    
    std::cout << "Generating OBJ mesh with " << vertices.size() << " vertices and "
      << quads.size() << " quads" << std::endl;
    
    // write vertices
    for(auto const & v : vertices) {
        file << "v " << v.x << ' ' << v.y << ' ' << v.z << '\n';
    }
    
    // write quad indices
    for(auto const & q : quads) {
        file << "f " << (q.i0+1) << ' ' << (q.i1+1) << ' ' << (q.i2+1) << ' ' << (q.i3+1) << '\n';
    }
    
    file.close();
}

//------------------------------------------------------------------------------

// Calculate the normal vector of a triangle of vertices.  The result isn't normalised.
void DualMCExample::triangleNormal(int v0, int v1, int v2, double &xn, double &yn, double &zn) const
{
	double x1 = vertices[v1].x - vertices[v0].x;
	double x2 = vertices[v2].x - vertices[v0].x;
	double y1 = vertices[v1].y - vertices[v0].y;
	double y2 = vertices[v2].y - vertices[v0].y;
	double z1 = vertices[v1].z - vertices[v0].z;
	double z2 = vertices[v2].z - vertices[v0].z;

	xn = y1*z2 - z1*y2;
	yn = z1*x2 - x1*z2;
	zn = y1*x2 - x1*y2;
}

// Write the quads out as pairs of triangles in an ASCII STL file
void DualMCExample::writeSTL(std::string const & fileName) const {
    std::cout << "Writing STL file " << fileName << std::endl;
    // check if we actually have an ISO surface
    if(vertices.size () == 0 || quads.size() == 0) {
        std::cout << "No ISO surface generated. Skipping STL generation." << std::endl;
        return;
    }
    
    // open output file
    std::ofstream file(fileName);
    if(!file) {
        std::cout << "Error opening output file" << std::endl;
        return;
    }
    
    std::cout << "Generating STL triangulation with " << vertices.size() << " vertices and "
      << 2*quads.size() << " triangles" << std::endl;

    file << "solid " << std::endl;
    
    double xn, yn, zn;

    // write quads as triangle pairs
    for(auto const & q : quads) 
    {
        triangleNormal(q.i0, q.i1, q.i2, xn, yn, zn);
        file << "facet normal " << xn << ' ' << yn << ' ' << zn << std::endl;
        file << " outer loop" << std::endl;
        file << "  vertex " << vertices[q.i0].x << ' ' << vertices[q.i0].y << ' ' << vertices[q.i0].z << std::endl;
        file << "  vertex " << vertices[q.i1].x << ' ' << vertices[q.i1].y << ' ' << vertices[q.i1].z << std::endl;
        file << "  vertex " << vertices[q.i2].x << ' ' << vertices[q.i2].y << ' ' << vertices[q.i2].z << std::endl;
        file << " endloop" << std::endl;
        file << "endfacet" << std::endl;

	triangleNormal(q.i0, q.i2, q.i3, xn, yn, zn);
        file << "facet normal " << xn << ' ' << yn << ' ' << zn << std::endl;
        file << " outer loop" << std::endl;
        file << "  vertex " << vertices[q.i0].x << ' ' << vertices[q.i0].y << ' ' << vertices[q.i0].z << std::endl;
        file << "  vertex " << vertices[q.i2].x << ' ' << vertices[q.i2].y << ' ' << vertices[q.i2].z << std::endl;
        file << "  vertex " << vertices[q.i3].x << ' ' << vertices[q.i3].y << ' ' << vertices[q.i3].z << std::endl;
        file << " endloop" << std::endl;
        file << "endfacet" << std::endl;
    }

    file << "endsolid " << std::endl;
    
    file.close();
}

//------------------------------------------------------------------------------

DualMCExample::RadialGaussian::RadialGaussian(
    float cX,
    float cY,
    float cZ,
    float variance
    ) : cX(cX), cY(cY), cZ(cZ) {
        float constexpr TWO_PI = 6.283185307179586f;
        normalization = 1.0f/sqrt(TWO_PI * variance);
        falloff = -0.5f / variance;
    }

//------------------------------------------------------------------------------

float DualMCExample::RadialGaussian::eval(float x, float y, float z) const {
    // compute squared input point distance to gauss center
    float const dx = x - cX;
    float const dy = y - cY;
    float const dz = z - cZ;
    float const dSquared = dx * dx + dy * dy + dz * dz;
    // compute gauss 
    return normalization * exp(falloff * dSquared);
}
