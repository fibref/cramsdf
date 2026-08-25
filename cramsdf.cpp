#include <numbers>
#include <algorithm>
#include <iostream>

#include "cramsdf.h"

using namespace cramsdf;

int solve_quadratic(float a, float b, float c, float x[2]) {
    float h = -b / 2;
    float discriminant = h * h - a * c;
    if (discriminant < 0) {
        return 0; // No real roots
    } else {
        float disc_sqrt = sqrt(discriminant);
        x[0] = (h - disc_sqrt) / a;
        x[1] = (h + disc_sqrt) / a;
        return 2; // Two real roots
    }
}

int solve_cubic(float a, float b, float c, float d, float x[3]) {
    if (a == 0)
        return solve_quadratic(b, c, d, x);
    b /= a; c /= a; d /= a;
    if (fabs(b) > 1e6)
        return solve_quadratic(b, c, d, x);

    // For equation x^3 + px + q = 0, let m = -p / 3, n = -q / 2.
    float bp2 = b * b;
    float m = 1.0f / 9.0f * (bp2 - 3 * c);
    float n = 1.0f / 54.0f * (9 * b * c - 27 * d - 2 * bp2 * b);
    float m3 = m * m * m;
    float discriminant = n * n - m3;
    float offset = -b / 3;

    if (discriminant > 0) {
        // float u = cbrt(n + sqrt(discriminant));
        // Somehow cbrt is a lot slower!
        float temp = n + sqrt(discriminant);
        float u = (temp < 0 ? -1 : 1) * pow(fabs(temp), 1.0f / 3.0f);
        float v = m / u;
        x[0] = u + v + offset;
        return 1;
    } else {
        float r = sqrt(m3);
        float temp = n / r;
        if (temp < -1) temp = -1;
        if (temp > 1) temp = 1;

        float theta = std::acos(temp) / 3;
        float amp = 2 * sqrt(m);

        constexpr float pi = std::numbers::pi_v<float>;
        x[0] = amp * std::cos(theta) + offset;
        x[1] = amp * std::cos(theta + 2 * pi / 3) + offset;
        x[2] = amp * std::cos(theta - 2 * pi / 3) + offset;
        return 3;
    }
}

std::pair<float, float> Segment::nearest_on_line(const Vec2& point) const {
    Vec2 to = p3 - p0;
    Vec2 ap = point - p0;
    Vec2 bp = point - p3;
    float t = dot(ap, to) / dot(to, to);
    if (t < 0) return { 0.0f, ap.len2() };
    if (t > 1) return { 1.0f, bp.len2() };

    Vec2 pq = to * t - ap;
    return { t, pq.len2() };
}

std::pair<float, float> Segment::nearest_on_conic(const Vec2& point) const {
    Vec2 delta1 = p1 - p0;
    Vec2 delta2 = p3 - p2 - delta1;

    Vec2 pa = p0 - point;
    Vec2 pc = p3 - point;
    float x[3];
    int num_roots = solve_cubic(
        dot(delta2, delta2),
        3 * dot(delta1, delta2),
        2 * dot(delta1, delta1) + dot(pa, delta2),
        dot(pa, delta1), x
    );

    float t_min = pa.len2() < pc.len2() ? 0.0f : 1.0f;
    float d2_min = t_min == 0.0f ? pa.len2() : pc.len2();

    for (int i = 0; i < num_roots; ++i) {
        if (x[i] < 0 || x[i] > 1) continue;
        Vec2 pq = delta2 * x[i] * x[i] + delta1 * 2 * x[i] + pa;
        if (pq.len2() < d2_min) {
            d2_min = pq.len2();
            t_min = x[i];
        }
    }
    return { t_min, d2_min };
}

std::pair<float, float> Segment::nearest_on_cubic(const Vec2& point) const {
    Vec2 pa = p0 - point;
    Vec2 pd = p3 - point;
    Vec2 delta1 = p1 - p0;
    Vec2 delta2 = p2 - p1 - delta1;
    Vec2 delta3 = (p3 - p2) - (p2 - p1) - delta2;

    float c0 = dot(pa, delta1) * 3;
    float c1 = dot(pa, delta2) * 6 + dot(delta1, delta1) * 9;
    float c2 = dot(pa, delta3) * 3 + dot(delta1, delta2) * 27;
    float c3 = dot(delta1, delta3) * 12 + dot(delta2, delta2) * 18;
    float c4 = dot(delta2, delta3) * 15;
    float c5 = dot(delta3, delta3) * 3;

    float t_min = pa.len2() < pd.len2() ? 0.0f : 1.0f;
    float d2_min = t_min == 0.0f ? pa.len2() : pd.len2();

    float t = 0.0f;
    for (int i = 0; i <= 4; ++i) {
        t = 1.0f / 4 * i;
        for (int j = 0; j <= 4; ++j) {
            float v = ((((c5 * t + c4) * t + c3) * t + c2) * t + c1) * t + c0;
            float dv = (((5 * c5 * t + 4 * c4) * t + 3 * c3) * t + 2 * c2) * t + c1;
            t = t - v / dv;
        }
        if (t < 0 || t > 1) continue;
        Vec2 pq = ((delta3 * t + delta2 * 3) * t + delta1 * 3) * t + pa;
        if (pq.len2() < d2_min) {
            d2_min = pq.len2();
            t_min = t;
        }
    }
    return { t_min, d2_min };
}

Vec2 Segment::at(float t, SegmentType type) const {
    switch (type) {
    case Line: return p0 + (p3 - p0) * t;
    case Conic: return (p3 - p2 - p1 + p0) * t * t + (p1 - p0) * 2 * t + p0;
    case Cubic: return (p3 - p2 * 3 + p1 * 3 - p0) * t * t * t + (p2 - p1 * 2 + p0) * 3 * t * t + (p1 - p0) * 3 * t + p0;
    default: return { 0.0f, 0.0f }; // Should never reach here.
    }
}
Vec2 Segment::tangent(float t, SegmentType type) const {
    switch (type) {
    case Line: return p3 - p0;
    case Conic: return (p3 - p2 - p1 + p0) * 2 * t + (p1 - p0) * 2;
    case Cubic: return (p3 - p2 * 3 + p1 * 3 - p0) * 3 * t * t + (p2 - p1 * 2 + p0) * 6 * t + (p1 - p0) * 3;
    default: return { 0.0f, 0.0f }; // Should never reach here.
    }
}

float Shape::signed_dist(const Vec2& point) const {
    float min_dist2 = std::numeric_limits<float>::max();
    float min_t = 0.0f;
    unsigned min_index = 0;

    for (unsigned i = 0; i < segments.size(); ++i) {
        auto [t, dist2] = segments[i].dist2(point, types[i]);
        if (dist2 < min_dist2) {
            min_dist2 = dist2;
            min_t = t;
            min_index = i;
        }
    }
    if (min_t > 0.0f && min_t < 1.0f) {
        Vec2 pq = segments[min_index].at(min_t, types[min_index]) - point;
        Vec2 tg = segments[min_index].tangent(min_t, types[min_index]);
        bool outside = (cross(pq, tg) > 0) ^ winding_order;
        return outside ? sqrt(min_dist2) : -sqrt(min_dist2);
    } else { // Nearest point is at a segment's endpoint.
        auto [a, b] = min_t == 0.0f ?
            std::pair(segments[prev(min_index)].exit(), segments[min_index].entry()) :
            std::pair(segments[min_index].exit(), segments[next(min_index)].entry());
        if (-dot(a, b) < std::cos(3) * a.len() * b.len()) {
            Vec2 q = min_t == 0.0f ? segments[min_index].p0 : segments[min_index].p3;
            Vec2 pq = q - point;
            Vec2 tg = a + b;
            bool outside = (cross(pq, tg) > 0) ^ winding_order;
            return outside ? sqrt(min_dist2) : -sqrt(min_dist2);
        }
        bool is_convex = (cross(a, b) < 0) ^ winding_order;
        return is_convex ? sqrt(min_dist2) : -sqrt(min_dist2);
    }
}

std::pair<Vec2, Vec2> Shape::cbounds() const {
    float x_min = std::numeric_limits<float>::max();
    float x_max = std::numeric_limits<float>::lowest();
    float y_min = std::numeric_limits<float>::max();
    float y_max = std::numeric_limits<float>::lowest();

    for (const Vec2* p = &segments[0].p0; p <= &segments.back().p3; ++p) {
        x_min = std::min(x_min, p->x);
        x_max = std::max(x_max, p->x);
        y_min = std::min(y_min, p->y);
        y_max = std::max(y_max, p->y);
    }
    return { Vec2(x_min, y_min), Vec2(x_max, y_max) };
}

std::pair<unsigned, unsigned> Shape::buffer_info(unsigned spread) const {
    auto [bl, tr] = cbounds();
    unsigned width = static_cast<unsigned>(tr.x - bl.x) + 2 * spread;
    unsigned height = static_cast<unsigned>(tr.y - bl.y) + 2 * spread;
    return { width, height };
}

void Shape::sdf(unsigned spread, unsigned char* buffer, unsigned width, unsigned height) const {
    auto [bl, tr] = cbounds();
    float px_len = (tr.x - bl.x) / static_cast<float>(width - 2 * spread);

    Vec2 pen = Vec2(bl.x - spread * px_len, bl.y - spread * px_len) + Vec2(0.5f * px_len, 0.5f * px_len);
    for (unsigned i = 0; i < height; ++i) {
        pen.x = bl.x - spread * px_len;
        for (unsigned j = 0; j < width; ++j) {
            float d = signed_dist(pen);
            int value = std::clamp(static_cast<int>(d * 127.5f / (spread * px_len) + 128.0f), 0, 255);
            *buffer = static_cast<unsigned char>(value); ++buffer;
            pen.x += px_len;
        }
        pen.y += px_len;
    }
}

struct Angle {
    Vec2 origin;
    Vec2 entry;
    Vec2 exit;
    bool is_convex(WindingOrder winding_order) const {
        return (cross(entry, exit) < 0) ^ winding_order;
    }
};

// op: The vector from the angle's origin to the point being evaluated.
std::pair<float, float> resolve_angle(const Angle& angle, const Vec2& op, WindingOrder winding_order) {
    auto dist = [](const Vec2& base, const Vec2& v) {
        float product = dot(base, v);
        float dist2 = v.len2() - product * product / base.len2();
        return std::sqrt(dist2);
        };

    float d0 = dist(angle.entry, op);
    float d1 = dist(angle.exit, op);

    int section = ((cross(op, angle.entry) > 0) ^ winding_order) << 1 | (cross(op, angle.exit) > 0) ^ winding_order;
    switch (section) {
    case 0: // Diagonal
        return { d0, d1 };
        break;
    case 1: // Entry(convex), Exit(concave)
        return { d0, -d1 };
        break;
    case 2: // Exit(convex), Entry(concave)
        return { -d0, d1 };
        break;
    case 3: // Inside
        return { -d0, -d1 };
        break;
    }
}

void Shape::cmsdf(unsigned spread, float threshold, unsigned char* buffer, unsigned width, unsigned height) const {
    // Collect all angles that are sharper than the threshold angle.
    std::vector<Angle> angles;
    for (unsigned i = 0; i < segments.size(); ++i) {
        Vec2 a = segments[i].exit();
        Vec2 b = segments[next(i)].entry();
        if (-dot(a, b) > std::cos(threshold) * a.len() * b.len())
            angles.push_back({ segments[i].p3, a, b });
    }
    std::cout << "Found " << angles.size() << " sharp angles." << std::endl;

    auto [bl, tr] = cbounds();
    float px_len = (tr.x - bl.x) / static_cast<float>(width - 2 * spread);

    // Radius for corner calculation
    const float corner_radius = 1.3f * px_len;
    // Radius for storing LSB flag
    const float lsb_radius = 2.5f * px_len;

    Vec2 pen = Vec2(bl.x - spread * px_len, bl.y - spread * px_len) + Vec2(0.5f * px_len, 0.5f * px_len);
    for (unsigned i = 0; i < height; ++i) {
        pen.x = bl.x - spread * px_len;
        for (unsigned j = 0; j < width; ++j) {
            // For each pixel, check if it is close to an angle.
            for (const auto& angle : angles) {
                Vec2 op = pen - angle.origin;
                float r = op.len();
                if (r < lsb_radius) {
                    if (r < corner_radius) {
                        auto [d0, d1] = resolve_angle(angle, op, winding_order);
                        int v0 = std::clamp(static_cast<int>(d0 * 127.5f / (spread * px_len) + 128.0f), 0, 255);
                        int v1 = std::clamp(static_cast<int>(d1 * 127.5f / (spread * px_len) + 128.0f), 0, 255);
                        if (angle.is_convex(winding_order))
                            v0 &= ~1;
                        else
                            v0 |= 1;

                        *buffer = static_cast<unsigned char>(v0); ++buffer;
                        *buffer = static_cast<unsigned char>(v1); ++buffer;
                    } else {
                        float d = signed_dist(pen);
                        int value = std::clamp(static_cast<int>(d * 127.5f / (spread * px_len) + 128.0f), 0, 255);
                        int v0 = value;
                        if (angle.is_convex(winding_order))
                            v0 &= ~1;
                        else
                            v0 |= 1;
                        *buffer = static_cast<unsigned char>(v0); ++buffer;
                        *buffer = static_cast<unsigned char>(value); ++buffer;
                    }
                    goto next_iteration;
                }
            }
            // If not close to an angle, compute the signed distance normally.
            {
                float d = signed_dist(pen);
                int value = std::clamp(static_cast<int>(d * 127.5f / (spread * px_len) + 128.0f), 0, 255);
                *buffer = static_cast<unsigned char>(value); ++buffer;
                *buffer = static_cast<unsigned char>(value); ++buffer;
            }
next_iteration:
            pen.x += px_len;
        }
        pen.y += px_len;
    }
}
