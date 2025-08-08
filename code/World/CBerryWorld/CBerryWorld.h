//  MABE is a product of The Hintze Lab @ MSU
//     for general research information:
//         hintzelab.msu.edu
//     for MABE documentation:
//         github.com/Hintzelab/MABE/wiki
//
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//     to view the full license, visit:
//         github.com/Hintzelab/MABE/wiki/License

#pragma once    // directive to insure that this .h file is only included one time

#include <World/AbstractWorld.h> // AbstractWorld defines all the basic function templates for worlds
#include <string>
#include <memory> // shared_ptr
#include <map>

#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>

// point2d defines a 2d vector with addtion, subtraction, dot/scalar product(*)
// and cross product
// also included are distance functions and functions which return the signed
// angle between 2 point2ds (relitive to 0,0)

class CPoint { // a point class, useful for 2d worlds
public:
    double x;
    double y;

    const double pi = atan(1.0) * 4;

    CPoint() {
        x = 0.0;
        y = 0.0;
    }
    CPoint(double _x, double _y) : x(_x), y(_y) {} // construct with x,y

    void set(double _x, double _y) {
        x = _x;
        y = _y;
    }

    // print x and y of this
    void show() const {
        std::cout << x << "," << y << "\n";
    }

    std::string as_string() const {
        return(std::to_string(x)+","+std::to_string(y));
    }

    // scalar/dot product of this and other
    CPoint operator=(const CPoint & other) {
        this->x = other.x;
        this->y = other.y;
        return *this;
    }

    // scalar/dot product of this and other
    double operator*(const CPoint & other) const {
        return x * other.x + y * other.y;
    }

    // vector * scalar
    CPoint scale(const double & mag) const {
        CPoint newVect(x * mag, y * mag);
        return newVect;
    }

    // cross product of this and other
    double cross_prod(const CPoint & other) const {
        return x * other.y - y * other.x;
    }

    // add this point and other
    CPoint operator+(const CPoint & other) const {
        CPoint newVect(x + other.x, y + other.y);
        return newVect;
    }

    // subtract other from this point
    CPoint operator-(const CPoint & other) const {
        CPoint newVect;
        newVect.x = x - other.x;
        newVect.y = y - other.y;
        return newVect;
    }

    // scalar/dot product of this and other
    bool operator==(const CPoint & other) const {
        if (x == other.x && y == other.y) {
            return true;
        }
        return false;
    }

    // length between this point and 0,0
    double dist() const {
        return sqrt(x * x + y * y);
    }

    // length between this point and other
    double dist(const CPoint & other) const {
        return (*this - other).dist();
    }

    // return angle in radians between this point and other relative to origin (0,0)
    double angle_between_radian(const CPoint & other) const {
        if (abs(x - other.x) < .000001 &&
            abs(y - other.y) < .000001) { // vectors are effecvily the same
            return (0);
        }
        if (abs(x / other.x - y / other.y) <
            .000001) { // vectors are effecvily parallel
            if (((x > 0) == (other.x > 0)) &&
                ((y > 0) == (other.y > 0))) { // and are pointing the same way
                return (0);
            }
        }
        return (cross_prod(other) < 0 ? 1 : -1) *
            acos((*this) * other / (dist() * other.dist()));
    }

    // return angle in degrees between this point and other relative to origin (0,0)
    double angle_between_deg(CPoint & other) const {
        return angle_between_radian(other) / pi * 180;
    }

    // check if this point is in a convex polygon defined by a list of points
    bool InConvexPoly(const std::vector<CPoint>& polygon) const {
        int side = 0;
        int i, i2;
        for (i = 0, i2 = polygon.size() - 1; i < polygon.size(); i2 = i++) {
            //If point is in the polygon
            if (polygon[i].x == this->x && polygon[i].y == this->y) {
                return true;
            }

            //Form a segment between the i'th point
            double x1 = polygon[i].x;
            double y1 = polygon[i].y;

            double x2 = polygon[i2].x;
            double y2 = polygon[i2].y;

            //Compute the cross product
            double d = (x - x1) * (y2 - y1) - (y - y1) * (x2 - x1);

            //std::cout << x1 << "," << y1 << "    " << x2 << "," << y2 << "   test: " << x << "," << y << "     " << d << std::endl;
            if (side == 0) {
                side = (d > 0) ? 1 : -1;
            }
            else {
                if (d * side < 0) {
                    return false;
                }
            }

        }
        //If no change in direction, then on same side of all segments, and thus inside
        return true;
    }

    // check if this point is in a convex polygon defined by a list of points
    bool InConvex4Poly(const std::array<CPoint,4>& polygon) const {
        int side = 0;
        int i, i2;
        for (i = 0, i2 = 3; i < 4; i2 = i++) {
            //If point is in the polygon
            if (polygon[i].x == this->x && polygon[i].y == this->y) {
                return true;
            }

            //Form a segment between the i'th point
            double x1 = polygon[i].x;
            double y1 = polygon[i].y;

            double x2 = polygon[i2].x;
            double y2 = polygon[i2].y;

            //Compute the cross product
            double d = (x - x1) * (y2 - y1) - (y - y1) * (x2 - x1);

            //std::cout << x1 << "," << y1 << "    " << x2 << "," << y2 << "   test: " << x << "," << y << "     " << d << std::endl;
            if (side == 0) {
                side = (d > 0) ? 1 : -1;
            }
            else {
                if (d * side < 0) {
                    return false;
                }
            }

        }
        //If no change in direction, then on same side of all segments, and thus inside
        return true;
    }

    // check if this point is in a box aligned with the x and y
    bool InBox(const std::pair<CPoint, CPoint> & corners) const {
        return (corners.first.x - x > 0)
            && (x - corners.second.x > 0)
            && (corners.first.y - y > 0)
            && (y - corners.second.y > 0);
    }

    // check if this point is in a circle
    bool InCircle(CPoint center, double radius) const {
        return (dist(center) < radius);
    }

    //Rotate a p counterclockwise around 0,0 by a angle a (indegrees)
    CPoint rotate(const double & a) const {
        double qx = cos(a / (180 / pi)) * this->x - sin(a / (180 / pi)) * this->y;
        double qy = sin(a / (180 / pi)) * this->x + cos(a / (180 / pi)) * this->y;
        return CPoint(qx,qy);
    }

    CPoint round() const {
        return CPoint(int(x), int(y));
    }

};


// Vector2d wraps a vector<T> and provides x,y style access
// no error checking is provided for out of range errors
// internally this class uses R(ow) and C(olumn) (i.e. how the data is stored in
// the data vector)
// the user sees x,y where x = column, y = row

template <typename T> class CGrid {
    std::vector<T> data;
    int R, C;

    // get index into data vector for a given x,y
    inline int getIndex(int r, int c) const { return (r * C) + c; }

public:
    CGrid() {
        R = 0;
        C = 0;
    }
    // construct a vector of size x * y
    CGrid(int x, int y) : R(y), C(x) { data.resize(R * C); }

    CGrid(int x, int y, T value) : R(y), C(x) { data.resize(R * C, value); }

    void reset(int x, int y) {
        R = y;
        C = x;
        data.clear();
        data.resize(R * C);
    }

    void reset(int x, int y, T value) {
        R = y;
        C = x;
        data.clear();
        data.resize(R * C, value);
    }

    // overwrite this classes data (vector<T>) with data coppied from newData
    void assign(std::vector<T> newData) {
        if ((int)newData.size() != R * C) {
            std::cout << "  ERROR :: in Vector2d::assign() vector provided does not "
                "fit. provided vector is size "
                << newData.size() << " but Rows(" << R << ") * Columns(" << C
                << ") == " << R * C << ". Exitting." << std::endl;
            exit(1);
        }
        data = newData;
    }

    // provides access to value x,y can be l-value or r-value (i.e. used for
    // lookup of assignment)
    T& operator()(int x, int y) { 
        return data[getIndex(y, x)];
    }

    T& operator()(double x, double y) const {
        return data[getIndex((int)(y), (int)(x))];
    }

    T& operator()(std::pair<int, int> loc) const {
        return data[getIndex(loc.second, loc.first)];
    }

    T& operator()(std::pair<double, double> loc) const {
        return data[getIndex((int)(loc.second), (int)(loc.first))];
    }

    T& operator()(CPoint loc) { 
        return data[getIndex((int)loc.y, (int)loc.x)];
    }

    // show the contents of this Vector2d with index values, and x,y values
    void show() const {
        for (int r = 0; r < R; r++) {
            for (int c = 0; c < C; c++) {
                std::cout << getIndex(r, c) << " : " << c << "," << r << " : "
                    << data[getIndex(r, c)] << "\n";
            }
        }
    }

    // show the contents of this Vector2d in a grid
    void showGrid(int precision = -1) const {
        if (precision < 0) {
            for (int r = 0; r < R; r++) {
                for (int c = 0; c < C; c++) {
                    std::cout << data[getIndex(r, c)] << " ";
                }
                std::cout << "\n";
            }
        }
        else {
            for (int r = 0; r < R; r++) {
                for (int c = 0; c < C; c++) {
                    if (data[getIndex(r, c)] == 0) {
                        std::cout << std::setfill(' ') << std::setw((precision * 2) + 2)
                            << " ";
                    }
                    else {
                        std::cout << std::setfill(' ') << std::setw((precision * 2) + 1)
                            << std::fixed << std::setprecision(precision)
                            << data[getIndex(r, c)] << " ";
                    }
                }
                std::cout << "\n";
            }
        }
    }
    int x() { return C; }

    int y() { return R; }
};


///////////////////////////////////////////////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////

class CBerryWorld : public AbstractWorld {


public:
    // parameters for group and brain namespaces
    static std::shared_ptr<ParameterLink<int>> evaluationsPerGenerationPL;
    
    // a local variable used for faster access to the ParameterLink value
    int evaluationsPerGeneration;
    
    std::string groupName = "root::";
    std::string brainName = "root::";
    
    int worldX = 25;
    int worldY = 25;

    // a viewWindow is basiclly an arc (simplifyed to a 4 sided polygon) that is used with the critter vision
    std::array<CPoint, 4> makeViewWindow(const double& near, const double& far, const double& width, const double& offset) {
        // near;    // how far from agent does this window start
        // far;     // how far from agent does this window end
        // width;   // in degrees
        // offset;  // degrees off center to rotate Window
        std::array<CPoint, 4> retunList;
        CPoint nearPoint = CPoint(0, near);
        CPoint farPoint = CPoint(0, far);
        retunList[0] = nearPoint.rotate(offset + (width / 2.0));
        retunList[1] = farPoint.rotate(offset + (width / 2.0));
        retunList[2] = farPoint.rotate(offset + (width / -2.0));
        retunList[3] = nearPoint.rotate(offset + (width / -2.0));
        return(retunList);
    };

    std::array<CPoint, 4> rotateViewWindow(const std::array<CPoint, 4>& viewWindow, const double& angle, const CPoint loc) {
        std::array<CPoint, 4> retunList;
        retunList[0] = viewWindow[0].rotate(angle) + loc;
        retunList[1] = viewWindow[1].rotate(angle) + loc;
        retunList[2] = viewWindow[2].rotate(angle) + loc;
        retunList[3] = viewWindow[3].rotate(angle) + loc;
        return(retunList);
    };

    class Plant {
    public:
        long long ID;
        int localIndex = -1;
        static long long nextID;
        CPoint loc = CPoint(0.0, 0.0);

        std::shared_ptr<Organism> org;

        long long parent = -1;
        int offCount = 0;
        double water = 1;
        double energy = 1; // resource update modulated by leaf coverage
        bool alive = true;

        double rootDist = 1;
        double rootWidth = 1;
        double leafCoverage = 1;
        double stalkHeight = 0;

        Plant() {
            ID = nextID++;
        }
    };

    class Critter {
    public:
        long long ID;
        int localIndex = -1;
        static long long nextID;
        CPoint loc = CPoint(0.0,0.0);
        double facing = Random::getDouble(360.0);
        std::shared_ptr<Organism> org;

        long long parent = -1;
        int offCount = 0;
        double energy = 1;
        bool alive = true;

        double moveCnt = 0;
        double turnCnt = 0;
        double eatCnt = 0;

        Critter() {
            ID = nextID++;
        }
    };

    class subSector {
        CPoint loc; // x,y of this Sector
        std::pair<CPoint, CPoint> corners;
        double resourceLevel;
        double waterLevel;
    };

    class Sector {
    public:
        CPoint loc; // x,y of this Sector
        std::pair<CPoint, CPoint> corners; // min and max corners of this sector
        CGrid<subSector> subSectors; // low resolution water and resource map, for plants
        std::vector<std::shared_ptr<Plant>> plants;
        std::vector< std::shared_ptr<Critter>> critters;
    };

    //CGrid<Sector> world;
    //std::vector<std::shared_ptr<Plant>> plants;
    //std::vector< std::shared_ptr<Critter>> critters;


    //class sector:

    CBerryWorld(std::shared_ptr<ParametersTable> PT);
    virtual ~CBerryWorld() = default;

    virtual auto evaluate(std::map<std::string, std::shared_ptr<Group>>& /*groups*/, int /*analyze*/, int /*visualize*/, int /*debug*/) -> void override;

    virtual auto requiredGroups() ->std::unordered_map<std::string, std::unordered_set<std::string>> override;

};

