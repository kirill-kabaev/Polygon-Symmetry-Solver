/**
 * @file symmetry_test.cpp
 * @brief Verified solver for convex polygon symmetry detection.
 * @author Kabaev Kirill
 * @note Task assigned by Song Xindi (Shenzhen Honeycomb Tech)
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <string>

/**
 * @brief Mathematical constant Pi.
 */
constexpr double PI = 3.14159265358979323846;

/**
 * @brief Numerical tolerance threshold for geometric floating-point comparisons.
 */
constexpr double EPS = 1e-5;

/**
 * @struct Point
 * @brief Represents a point in 2D Cartesian coordinates.
 */
struct Point {
    double x; ///< Abscissa coordinate.
    double y; ///< Ordinate coordinate.

    /**
     * @brief Default constructor initializing point at the origin (0.0, 0.0).
     */
    Point() : x(0.0), y(0.0) {}

    /**
     * @brief Parameterized constructor.
     * @param _x Coordinate on X-axis.
     * @param _y Coordinate on Y-axis.
     */
    Point(double _x, double _y) : x(_x), y(_y) {}

    /**
     * @brief Equality comparison operator within tolerance EPS.
     * @param other Target point to compare against.
     * @return True if both coordinates differ by less than EPS, false otherwise.
     */
    bool operator==(const Point& other) const {
        return std::abs(x - other.x) < EPS && std::abs(y - other.y) < EPS;
    }
};

/**
 * @brief Computes Euclidean distance between two points in 2D space.
 * @param a First point.
 * @param b Second point.
 * @return Scalar Euclidean distance between point a and point b.
 */
double distance(Point a, Point b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return std::sqrt(dx * dx + dy * dy);
}

/**
 * @brief Computes the arithmetic centroid (center of mass) of polygon vertices.
 * @param pts Vector of polygon vertices.
 * @return Point representing the arithmetic mean of all coordinates.
 */
Point compute_centroid(const std::vector<Point>& pts) {
    double cx = 0.0;
    double cy = 0.0;
    double n = static_cast<double>(pts.size());
    for (const auto& p : pts) {
        cx += p.x;
        cy += p.y;
    }
    return Point(cx / n, cy / n);
}

/**
 * @brief Reflects a point across an infinite line defined by two distinct points.
 * @param p Point to be reflected.
 * @param p1 First point defining the reflection axis line.
 * @param p2 Second point defining the reflection axis line.
 * @return Reflected point coordinates in 2D space.
 */
Point reflect_point(Point p, Point p1, Point p2) {
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    double norm_sq = dx * dx + dy * dy;
    if (norm_sq < 1e-12) return p;

    double vx = p.x - p1.x;
    double vy = p.y - p1.y;
    double t = (vx * dx + vy * dy) / norm_sq;
    Point proj(p1.x + t * dx, p1.y + t * dy);
    return Point(2.0 * proj.x - p.x, 2.0 * proj.y - p.y);
}

/**
 * @brief Checks whether a line preserves the polygon under planar reflection.
 * @param polygon Ordered cyclic list of polygon vertices.
 * @param p1 First point on the candidate axis of symmetry.
 * @param p2 Second point on the candidate axis of symmetry.
 * @return True if reflected vertices match original vertices with reversed cyclic orientation.
 */
bool is_symmetry_axis(const std::vector<Point>& polygon, Point p1, Point p2) {
    const int n = static_cast<int>(polygon.size());
    std::vector<Point> reflected(n);
    for (int i = 0; i < n; ++i) {
        reflected[i] = reflect_point(polygon[i], p1, p2);
    }

    std::vector<int> mapped_indices(n, -1);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (distance(reflected[i], polygon[j]) < EPS) {
                mapped_indices[i] = j;
                break;
            }
        }
        if (mapped_indices[i] == -1) return false;
    }

    std::vector<bool> seen(n, false);
    for (int idx : mapped_indices) {
        if (idx < 0 || idx >= n || seen[idx]) return false;
        seen[idx] = true;
    }

    int start = mapped_indices[0];
    for (int i = 0; i < n; ++i) {
        int expected = (start - i + n) % n;
        if (mapped_indices[i] != expected) return false;
    }

    return true;
}

/**
 * @brief Checks whether point p lies on the infinite line passing through c and f.
 * @param p Evaluated point.
 * @param c First reference point on the line (centroid).
 * @param f Second reference point on the line (feature).
 * @return True if point p lies on line (c, f) within tolerance.
 */
bool is_on_line(Point p, Point c, Point f) {
    double cross = (p.x - c.x) * (f.y - c.y) - (p.y - c.y) * (f.x - c.x);
    return std::abs(cross) < 1e-6;
}

/**
 * @brief Searches for all unique reflectional symmetry axes of a convex polygon.
 * @param polygon Ordered cyclic list of polygon vertices.
 * @return Vector of pairs of boundary points defining each valid symmetry axis.
 */
std::vector<std::pair<Point, Point>> find_symmetry_axes(const std::vector<Point>& polygon) {
    const int n = static_cast<int>(polygon.size());
    if (n < 3) return {};

    Point centroid = compute_centroid(polygon);

    // Feature points: vertices and edge midpoints
    std::vector<Point> features;
    features.reserve(2 * n);
    for (int i = 0; i < n; ++i) {
        features.push_back(polygon[i]);
    }
    for (int i = 0; i < n; ++i) {
        features.push_back(Point(
            (polygon[i].x + polygon[(i + 1) % n].x) / 2.0,
            (polygon[i].y + polygon[(i + 1) % n].y) / 2.0
        ));
    }

    std::vector<std::pair<Point, Point>> valid_axes;
    std::vector<double> seen_angles;

    for (const auto& f : features) {
        if (distance(centroid, f) < 1e-6) continue;

        double dx = f.x - centroid.x;
        double dy = f.y - centroid.y;
        double ang = std::atan2(dy, dx);
        if (ang < 0.0) ang += PI;
        if (ang >= PI) ang -= PI;

        bool duplicate = false;
        for (double prev : seen_angles) {
            double diff = std::abs(ang - prev);
            if (diff < 1e-5 || std::abs(diff - PI) < 1e-5) {
                duplicate = true;
                break;
            }
        }
        if (duplicate) continue;

        if (is_symmetry_axis(polygon, centroid, f)) {
            seen_angles.push_back(ang);

            // Locate boundary points on opposite sides of centroid along the axis
            Point p_pos = f;
            Point p_neg = f;
            double max_pos_dist = -1.0;
            double max_neg_dist = -1.0;

            for (const auto& feat : features) {
                if (is_on_line(feat, centroid, f)) {
                    double dot = (feat.x - centroid.x) * dx + (feat.y - centroid.y) * dy;
                    double d = distance(feat, centroid);

                    if (dot > 1e-6 && d > max_pos_dist) {
                        max_pos_dist = d;
                        p_pos = feat;
                    } else if (dot < -1e-6 && d > max_neg_dist) {
                        max_neg_dist = d;
                        p_neg = feat;
                    }
                }
            }

            if (max_pos_dist > 0 && max_neg_dist > 0) {
                valid_axes.push_back({p_neg, p_pos});
            } else {
                valid_axes.push_back({centroid, f});
            }
        }
    }

    return valid_axes;
}

/**
 * @brief Prints a coordinate value omitting redundant decimal zeros for integers.
 * @param val Double value to print.
 */
void print_coord(double val) {
    if (std::abs(val) < 1e-7) val = 0.0;
    if (std::abs(val - std::round(val)) < 1e-6) {
        std::cout << static_cast<long long>(std::round(val));
    } else {
        std::cout << val;
    }
}

/**
 * @brief Outputs the detected symmetry axes or reports non-symmetric status.
 * @param axes Vector of pairs of boundary points defining symmetry axes.
 */
void print_axes(const std::vector<std::pair<Point, Point>>& axes) {
    if (axes.empty()) {
        std::cout << "non-symmetric\n";
    } else {
        for (const auto& ax : axes) {
            print_coord(ax.first.x); std::cout << " ";
            print_coord(ax.first.y); std::cout << " - ";
            print_coord(ax.second.x); std::cout << " ";
            print_coord(ax.second.y); std::cout << "\n";
        }
    }
}

/**
 * @brief Executes and prints verification output for a single polygon test case.
 * @param title Title of the test example.
 * @param poly Vector of polygon node coordinates.
 */
void run_example(const std::string& title, const std::vector<Point>& poly) {
    std::cout << title << "\n";
    std::cout << "Nodes:\n";
    for (const auto& p : poly) {
        print_coord(p.x);
        std::cout << " ";
        print_coord(p.y);
        std::cout << "\n";
    }

    std::cout << "Output:\n";
    auto axes = find_symmetry_axes(poly);
    print_axes(axes);
    std::cout << "\n";
}

/**
 * @brief Main entry point supporting both file processing (CLI) and demo tests.
 * @param argc Number of command-line arguments.
 * @param argv Command-line arguments array.
 * @return 0 on success, 1 on file reading or parsing error.
 */
int main(int argc, char* argv[]) {
    // 1. Run standard built-in examples first
    std::cout << "=== Running Built-in Verification Examples ===\n\n";

    // Example 1.
    std::vector<Point> ex1 = {
        {0, 0},
        {1, 0},
        {1, 1},
        {0, 1}
    };
    run_example("Example 1.", ex1);

    // Example 2.
    std::vector<Point> ex2 = {
        {0, 0},
        {2, 1},
        {0, 3},
        {-2, 1}
    };
    run_example("Example 2.", ex2);

    // Example 3.
    std::vector<Point> ex3 = {
        {0.1, 1.0},
        {-1.0, 0.0},
        {0.0, -1.0},
        {1.0, -0.5},
        {2.0, 1.0}
    };
    run_example("Example 3.", ex3);

    // 2. Process file from command line options if provided
    std::cout << "=== Command Line File Execution ===\n";
    if (argc < 2) {
        std::cout << "Input (command line options): no text file specified.\n";
        std::cout << "Usage: " << argv[0] << " <filename.txt>\n";
        return 0;
    }

    std::string filename = argv[1];
    std::cout << "Input (command line options): " << filename << "\n";

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file '" << filename << "'\n";
        return 1;
    }

    std::vector<Point> user_polygon;
    double x, y;
    while (file >> x >> y) {
        user_polygon.push_back(Point(x, y));
    }
    file.close();

    if (user_polygon.size() < 3) {
        std::cerr << "Error: Polygon in file must contain at least 3 vertices.\n";
        return 1;
    }

    std::cout << "Output:\n";
    auto file_axes = find_symmetry_axes(user_polygon);
    print_axes(file_axes);

    return 0;
}
