// Copyright (C) 2018, Dominik Wodniok
// This software may be modified and distributed under the terms
// of the BSD 3-Clause license.
// See the LICENSE.txt file for details.

#ifndef GENTABLES_H_INCLUDED
#define GENTABLES_H_INCLUDED

#include <cstdint>
#include <map>
#include <vector>

/// Application class for generating the dual marching cubes and manifold
/// dual marching cubes tables.
class GenerateTablesApp {
public:
    /// Run the application.
    void run();
    
private:
    // forward declarations
    class CoordinateAxis;
    class CubeConfiguration;
    class CubeCornerCode;
    
    // copy from dualmc.h
    enum DMCEdgeCode {
        EDGE0 = 1,
        EDGE1 = 1 << 1,
        EDGE2 = 1 << 2,
        EDGE3 = 1 << 3,
        EDGE4 = 1 << 4,
        EDGE5 = 1 << 5,
        EDGE6 = 1 << 6,
        EDGE7 = 1 << 7,
        EDGE8 = 1 << 8,
        EDGE9 = 1 << 9,
        EDGE10 = 1 << 10,
        EDGE11 = 1 << 11,
        FORCE_32BIT = 0xffffffff
    };
    
    using ProblematicConfigsMap = std::map<uint8_t,uint32_t>;
    
private:
    // functions for generating the manifold dual marching cubes table
    template<class T> void registerConfigAxisRotations(T f, CubeConfiguration cubeConfiguration, CoordinateAxis const axis, ProblematicConfigsMap & problematicConfigsMap);
    void exploreConfigRotations(CubeConfiguration config, ProblematicConfigsMap & problematicConfigsMap);
    void generateManifoldTable(ProblematicConfigsMap & problematicConfigsMap);
    void writeManifoldTable(ProblematicConfigsMap const & problematicConfigsMap);
    
    // functions for generating the dual marching cubes table
    void generateDualMCTable(std::vector<uint32_t> & dualPointsList);
    void writeDualMCTable(std::vector<uint32_t> const & dualPointsList);
    
    /// Table for retrieving the edges adjacent to a cube corner.
    static uint32_t const cornerEdges[8][3];
};

/// Class for a cube configuration represented by the corners in/out classification.
class GenerateTablesApp::CubeConfiguration {
public:
    // in bit-mask for for each corner
    static constexpr uint8_t C0 = 1;
    static constexpr uint8_t C1 = 2;
    static constexpr uint8_t C2 = 4;
    static constexpr uint8_t C3 = 8;
    static constexpr uint8_t C4 = 16;
    static constexpr uint8_t C5 = 32;
    static constexpr uint8_t C6 = 64;
    static constexpr uint8_t C7 = 128;
    
    /// Initializing constructor.
    CubeConfiguration(uint8_t c) : config(c) {}
    
    /// Get the configuration value
    uint8_t get() const { return config; }
    
    /// Rotate cube configuration around x axis.
    void rotX();
    
    /// Rotate cube configuration around y axis.
    void rotY();
    
    /// Rotate cube configuration around z axis.
    void rotZ();
private:
    uint8_t config;
};

/// Class for representing a coordinate axis.
class GenerateTablesApp::CoordinateAxis {
public:
    /// Axis value enum
    enum class Value : uint32_t {
        NX = 0,
        PX = 1,
        NY = 2,
        PY = 3,
        NZ = 4,
        PZ = 5
    };
    
    /// Initializing constructor
    CoordinateAxis(Value v);
    
    /// Rotate axis around x axis.
    void rotX();
    
    /// Rotate axis around y axis.
    void rotY();
    
    /// Rotate axis around z axis.
    void rotZ();
    
    /// Get the stored axis value.
    Value get() const;
private:
    Value value;
};

/// Cuber corner stored as a Morton node.
class GenerateTablesApp::CubeCornerCode {
public:
    /// Initializing constructor.
    CubeCornerCode(uint8_t c) : code(c) {}
    
    /// Get the bit mask where the bit corresponding to the code value is set.
    uint32_t getMask() const { return 1u << code;}
    
    /// Get the Morton code.
    uint32_t getCode() const { return code;}
    
    /// Get the code of the neighbor in x direction.
    uint32_t nX() const { return code ^ 1;}
    
    /// Get the code of the neighbor in y direction.
    uint32_t nY() const { return code ^ 2;}
    
    /// Get the code of the neighbor in z direction.
    uint32_t nZ() const { return code ^ 4;}
private:
    /// Morton code of the corner
    uint32_t code;
};

#endif // GENTABLES_H_INCLUDED