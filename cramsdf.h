#pragma once

#include <utility>
#include <vector>
#include <cmath>

namespace cramsdf {

class Vec2 {
public:
    float x, y;

    Vec2 operator+(const Vec2& other) const {
        return { x + other.x, y + other.y };
    }
    Vec2 operator-(const Vec2& other) const {
        return { x - other.x, y - other.y };
    }
    Vec2 operator*(float scalar) const {
        return { x * scalar, y * scalar };
    }
    bool operator==(const Vec2& other) const = default;
    float len() const {
        return std::sqrt(x * x + y * y);
    }
    float len2() const {
        return x * x + y * y;
    }

};
inline float dot(const Vec2& a, const Vec2& b) {
    return a.x * b.x + a.y * b.y;
}
inline float cross(const Vec2& a, const Vec2& b) {
    return a.x * b.y - a.y * b.x;
}

enum SegmentType { Line, Conic, Cubic };

// The winding order of the outer contour.
enum WindingOrder { CW = 0, CCW = 1 };

class Segment {
public:
    Vec2 p0, p1, p2, p3;

public:
    Segment(const Vec2& start, const Vec2& end)
        : p0(start)
        , p1(end)
        , p2(start)
        , p3(end) {}

    Segment(const Vec2& start, const Vec2& c, const Vec2& end)
        : p0(start)
        , p1(c)
        , p2(c)
        , p3(end) {}

    Segment(const Vec2& start, const Vec2& c1, const Vec2& c2, const Vec2& end)
        : p0(start)
        , p1(c1)
        , p2(c2)
        , p3(end) {}

    Vec2 entry() const {
        return p1 - p0;
    }
    Vec2 exit() const {
        return p3 - p2;
    }
    std::pair<float, float> dist2(const Vec2& point, SegmentType type) const {
        switch (type) {
        case Line: return nearest_on_line(point);
        case Conic: return nearest_on_conic(point);
        case Cubic: return nearest_on_cubic(point);
        default: return { 0.0f, 0.0f }; // Should never reach here.
        }
    }
    Vec2 at(float t, SegmentType type) const;
    Vec2 tangent(float t, SegmentType type) const;

private:
    std::pair<float, float> nearest_on_line(const Vec2& point) const;
    std::pair<float, float> nearest_on_conic(const Vec2& point) const;
    std::pair<float, float> nearest_on_cubic(const Vec2& point) const;
};

class Shape {
public:
    std::vector<Segment> segments;
    std::vector<SegmentType> types;
    // Stores the indices of the segments that start a contour.
    std::vector<unsigned> contours;

    float signed_dist(const Vec2& point, WindingOrder winding_order) const;
    std::pair<Vec2, Vec2> cbounds() const;
    std::pair<unsigned, unsigned> buffer_info(unsigned spread) const;
    void sdf(unsigned spread, unsigned char* buffer, unsigned width, unsigned height) const;

private:
    unsigned next(unsigned index) const {
        for (auto it = contours.cbegin(); it != contours.cend(); ++it) {
            if (index < *it - 1) return index + 1;
            if (index == *it - 1) return *(it - 1);
        }
        return 1; // Should never reach here.
    }
    unsigned prev(unsigned index) const {
        for (auto it = contours.crbegin(); it != contours.crend(); ++it) {
            if (index > *it) return index - 1;
            if (index == *it) return *(it - 1) - 1;
        }
        return 1; // Should never reach here.
    }
};

}