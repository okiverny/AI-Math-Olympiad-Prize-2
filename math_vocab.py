'''
Collection of Math words
'''

math_words = []

## Algebra
# 1. Polynomials & Factorization – Factor theorem, Vieta’s formulas, symmetric sums.
algebra_polynomials_words = [
    "polynomial", "degree", "coefficient", "root", "factor", "theorem", "equation",
    "expression", "variable", "monomial", "binomial", "trinomial", "quadratic", 
    "cubic", "quartic", "quintic", "symmetric", "sum", "product", "identity", 
    "integer", "rational", "irrational", "real", "complex", "algebraic", "zero", 
    "constant", "term", "solution", "substitution", "discriminant", "remainder", 
    "divisibility", "synthetic", "division", "gcd", "lcm", "irreducible", "decomposition",
    "transform", "homogeneous", "nonlinear", "linear", "sequence", "series", "relation",
    "binomial", "expansion", "newton", "sums", "power", "degree", "leading", "term",
    "derivative", "differentiation", "integration", "substitution", "recurrence",
    "symmetry", "alternating", "alternant", "permutation", "combination", "factorial",
    "modular", "modulus", "prime", "composite", "congruence", "diophantine", "sequence",
    "finite", "infinite", "recursive", "closed", "roots", "radical", "simplification",
    "substitution", "transformation", "mapping", "bijection", "surjection", "injection",
    "determinant", "matrix", "trace", "characteristic", "invariant", "eigenvalue",
    "eigenvector", "canonical", "quadratic", "cubic", "logarithm", "exponent", 
    "multiplicity", "reduction", "complex number", "real number", "imaginary number",
    "Vieta", "Newton", "Euler", "Lagrange", "Gauss", "Bezout", "Fermat",
    "Cauchy", "Galois", "Descartes", "Taylor", "Maclaurin", "Bernoulli",
    "Chebyshev", "Kronecker", "Ruffini", "Abel", "Germain", "Laplace", "Sylvester"
]
math_words.extend(algebra_polynomials_words)

# 2. Equations & Inequalities – Quadratic equations, functional equations, AM-GM, Cauchy-Schwarz.
algebra_equations_inequalities_words = [
    # Quadratic Equations
    "quadratic", "roots", "discriminant", "parabola", "vertex",
    "symmetry", "axis", "focus", "directrix", "coefficients",
    "solution", "factorization", "Vieta", "sum", "product",
    "polynomial", "irreducible", "real", "complex", "imaginary",

    # Functional Equations
    "function", "mapping", "bijection", "involution", "recursion",
    "injective", "surjective", "bijective", "symmetry", "inverse",
    "composition", "substitution", "identity", "iteration", "transformation",
    "homogeneous", "nonlinear", "domain", "range", "codomain",

    # Inequalities (AM-GM, Cauchy-Schwarz, Jensen’s)
    "inequality", "AM-GM", "Cauchy-Schwarz", "Jensen", "Muirhead",
    "triangle", "arithmetic", "geometric", "harmonic", "convexity",
    "concave", "weighted", "optimization", "bound", "maximum",
    "minimum", "extremal", "absolute", "relative", "mean",

    # Algebraic Properties and Number Theory
    "expression", "identity", "symmetry", "modulus", "absolute",
    "rational", "irrational", "integer", "divisibility", "greatest",
    "least", "multiple", "prime", "composite", "coprime",
    "modulo", "exponent", "logarithm", "radical", "exponential",

    # Equations (General)
    "equation", "solution", "unknown", "system", "substitution",
    "elimination", "homogeneous", "inhomogeneous", "parametric", "determinant",
    "matrix", "vector", "linear", "nonlinear", "cubic",
    "quartic", "binomial", "logarithmic", "trigonometric", "exponential",

    # Names
    "Vieta", "Euler", "Lagrange", "Cauchy", "Schwarz",
    "Jensen", "Muirhead", "Chebyshev", "Bernoulli", "Gauss",
    "Newton", "Descartes", "Hilbert", "Laplace", "Abel",
    "Fermat", "Sylvester", "Noether", "Riemann", "Kronecker"
]
math_words.extend(algebra_equations_inequalities_words)

# 3. Sequences & Series – Arithmetic and geometric progressions, telescoping sums, recurrence relations.
algebra_sequences_words = [
    # Arithmetic and Geometric Progressions
    "sequence", "series", "term", "progression", "arithmetic",
    "geometric", "common", "difference", "ratio", "sum",
    "partial", "finite", "infinite", "first", "last",
    "index", "nth", "formula", "derivation", "general",

    # Telescoping Sums
    "telescoping", "collapsing", "cancellation", "simplification", "summation",
    "pairwise", "identity", "finite", "expansion", "reduction",
    "factorization", "rearrangement", "trivial", "substitution", "manipulation",
    "partition", "convergence", "divergence", "approximation", "alternating",

    # Recurrence Relations
    "recurrence", "relation", "recursive", "initial", "condition",
    "homogeneous", "nonhomogeneous", "solution", "induction", "characteristic",
    "roots", "polynomial", "linear", "nonlinear", "order",
    "closed-form", "iteration", "derivation", "proof", "explicit",

    # Summation Notation and Techniques
    "summation", "notation", "sigma", "index", "bounds",
    "upper", "lower", "finite", "infinite", "convergent",
    "divergent", "approximate", "estimation", "series", "expansion",
    "product", "factorial", "binomial", "harmonic", "Euler",

    # Special Sequences & Series
    "harmonic", "Fibonacci", "Lucas", "Catalan", "Bernoulli",
    "binomial", "triangular", "pentagonal", "power", "Taylor",
    "Maclaurin", "Newton", "Lagrange", "Chebyshev", "Legendre",
    "Dirichlet", "Stirling", "Eulerian", "hypergeometric", "generating",

    "Euler", "Fibonacci", "Lucas", "Bernoulli", "Taylor",
    "Maclaurin", "Newton", "Lagrange", "Chebyshev", "Legendre",
    "Dirichlet", "Stirling", "Vandermonde", "Cauchy", "Kronecker",
    "Ramanujan", "Gauss", "Laplace", "Pascal", "Catalan"
]
math_words.extend(algebra_sequences_words)

# 4. Number Theory – Modular arithmetic, prime numbers, divisibility rules, Diophantine equations.
algebra_sequences_words = [
    # Modular Arithmetic
    "modular", "modulo", "congruence", "residue", "inverse",
    "mod", "remainder", "division", "divisibility", "modulus",
    "equivalence", "order", "primitive", "power", "exponentiation",
    "cyclic", "subgroup", "inverse", "least", "greatest",

    # Prime Numbers
    "prime", "composite", "coprime", "relatively", "gcd",
    "lcm", "factor", "factorization", "sieve", "totient",
    "perfect", "square-free", "divisor", "multiplicative", "function",
    "primality", "test", "smallest", "largest", "density",

    # Divisibility Rules
    "divisibility", "division", "multiple", "remainder", "quotient",
    "divisor", "unit", "greatest", "common", "least",
    "sum", "difference", "product", "even", "odd",
    "binary", "decimal", "representation", "consecutive", "invariant",

    # Diophantine Equations
    "Diophantine", "integer", "solution", "equation", "linear",
    "Pell", "quadratic", "sum", "difference", "positive",
    "negative", "exponent", "radical", "fraction", "rational",
    "irrational", "integral", "square", "cube", "power",

    # Special Topics in Number Theory
    "Euler", "Fermat", "Wilson", "Chinese", "Remainder",
    "Legendre", "Jacobi", "order", "multiplicative", "Fibonacci",
    "reciprocal", "continued", "fraction", "Bezout", "conjecture",
    "theorem", "lemma", "corollary", "proof", "induction",

    "Euler", "Fermat", "Gauss", "Legendre", "Dirichlet",
    "Lagrange", "Ramanujan", "Carmichael", "Pell", "Pythagoras",
    "Mersenne", "Sophie", "Lucas", "Vandermonde", "Chebyshev",
    "Kronecker", "Pascal", "Newton", "Catalan", "Waring"
]
math_words.extend(algebra_sequences_words)

# 5. Combinatorial Algebra – Summation identities, generating functions, binomial coefficients.
algebra_combinatorics_words = [
    # Summation Identities
    "summation", "series", "sum", "sigma", "telescoping",
    "index", "bound", "partial", "finite", "infinite",
    "sequence", "recurrence", "formula", "closed", "expression",
    "expansion", "rearrangement", "manipulation", "inversion", "substitution",
    
    # Generating Functions
    "generating", "function", "power", "series", "polynomial",
    "coefficient", "recursion", "exponential", "ordinary", "formal",
    "convolution", "transform", "product", "growth", "roots",
    "factorial", "asymptotic", "partition", "convergence", "representation",

    # Binomial Coefficients
    "binomial", "coefficient", "choose", "combinatorial", "identity",
    "Newton", "expansion", "Pascal", "row", "triangle",
    "hockey-stick", "multinomial", "distribution", "formula", "summand",
    "indexing", "lattice", "weight", "catalan", "stirling",

    # Other Combinatorial Algebra Topics
    "counting", "arrangement", "bijection", "bijection", "bijection",
    "composition", "decomposition", "enumeration", "fibonacci", "harmonic",
    "alternating", "symmetry", "orthogonality", "approximation", "inequality",
    "divisibility", "modulo", "subset", "superset", "partition",
    "transformation", "functional", "coefficient-wise", "induction", "proof",
    "expansion", "manipulation", "relation", "generality", "specialization",
    "homogeneous", "heterogeneous", "equality", "inequality", "bijection",

    "Euler", "Newton", "Vandermonde", "Pascal", "Lagrange",
    "Fermat", "Stirling", "Catalan", "Bernoulli", "Ramanujan",
    "Chebyshev", "Dirichlet", "Kronecker", "Gauss", "Jacobi",
    "Riemann", "Taylor", "Wallis", "Cauchy", "Titchmarsh"
]
math_words.extend(algebra_combinatorics_words)

## Geometry
# 6. Plane Geometry – Triangles, circles, cyclic quadrilaterals, similarity, congruence.
plane_geometry_words = [
    "triangle", "circle", "quadrilateral", "cyclic", "similarity", "congruence", "angle", "perpendicular", "bisector", "median",
    "altitude", "orthocenter", "centroid", "circumcenter", "incenter", "radius", "diameter", "chord", "tangent", "secant",
    "arc", "sector", "segment", "inscribed", "circumscribed", "excircle", "incircle", "midpoint", "parallel", "intersection",
    "perpendicular bisector", "angle bisector", "Pythagorean theorem", "trigonometry", "sine", "cosine", "tangent", "cotangent",
    "secant line", "tangent line", "interior angle", "exterior angle", "sum of angles", "isosceles", "equilateral", "scalene",
    "right triangle", "acute", "obtuse", "altitude foot", "Euler line", "nine-point circle", "Brocard point", "Steiner line",
    "Feuerbach circle", "Simson line", "Menelaus theorem", "Ceva's theorem", "power of a point", "radical axis", "radical center",
    "homothety", "dilation", "reflection", "translation", "rotation", "symmetry", "harmonic division", "cross ratio", "perspective",
    "pole", "polar", "desargues theorem", "Pappus theorem", "Pascal's theorem", "Brianchon's theorem", "incircle-excircle theorem",
    "equal angles", "concurrent lines", "collinear points", "circumradius", "inradius", "semi-perimeter", "area", "heron's formula",
    "Brahmagupta's formula", "diagonal", "trapezoid", "parallelogram", "rhombus", "rectangle", "square", "kite",
    "cyclic quadrilateral theorem", "cyclicity", "opposite angles", "equal chords", "tangent-secant theorem", "butterfly theorem",
    "tangential quadrilateral", "bicentric quadrilateral", "perimeter", "Euler's theorem", "Apollonius theorem", "Thales theorem",
    "Simson theorem", "Miquel theorem", "Stewart's theorem", "Mordell theorem", "Archimedean property", "tangent circles",
    "nested circles", "concyclic points", "tangent-tangent theorem", "circle inversion", "pole-polar theorem", "tangent chord angle",
    "curvilinear triangle", "arc midpoint", "tangency point", "Ptolemy's theorem", "Newton's theorem", "inversive geometry",
    "constructible points", "golden ratio", "parallel postulate", "taxicab geometry", "Minkowski space", "affine geometry",
    "projective geometry", "hyperbolic geometry", "spherical geometry", "convex hull", "barycentric coordinates",
    "Carlyle circle", "Morley's theorem", "Steiner chain", "Johnson's theorem", "locus", "geometric mean", "Fermat's point",
    "Napoleon's theorem", "Van Aubel's theorem", "Lemoine's theorem", "Harmonic quadrilateral", "pentagon", "hexagon",
    "tessellation", "Hessian determinant", "Schläfli symbol", "Gergonne point", "Nagel point", "Mittenpunkt", "Malfatti circles",
    "Hagge circle", "extracenters"
]
math_words.extend(plane_geometry_words)

# 7. Coordinate Geometry – Lines, conics, transformations, equations of curves.
coordinate_geometry_words = [
    "line", "slope", "intercept", "point", "distance", 
    "midpoint", "perpendicular", "parallel", "segment", "vector", 
    "coordinates", "origin", "axes", "slope-intercept form", "point-slope form",
    "distance formula", "angle bisector", "dot product", "cross product", "vector addition", 
    "equation of a line", "line segment", "collinearity", "direction ratios", "normal vector",
    "orthogonal", "reflection", "rotation", "translation", "scaling", 
    "matrix", "determinant", "affine transformation", "homogeneous coordinates", "conic", 
    "parabola", "ellipse", "hyperbola", "circle", "focus", "directrix", "latus rectum", 
    "vertex", "axis of symmetry", "eccentricity", "standard form", "general form", 
    "polar coordinates", "conic section", "parabola equation", "ellipse equation", 
    "hyperbola equation", "circle equation", "circle center", "circle radius", "foci", 
    "asymptotes", "symmetry", "axis of symmetry", "tangent", "secant", "normal line", 
    "intersection", "angle between lines", "distance from a point to a line", 
    "equation of a circle", "circle intersection", "transformation matrix", 
    "linear transformation", "quadratic form", "inverse transformation", "projective geometry", 
    "coordinates transformation", "coordinate rotation", "homogeneous transformation", 
    "affine geometry", "convex hull", "polar form", "parametric equations", "curve", 
    "parameter", "line equation", "symmetric equations", "implicit equation", "explicit equation", 
    "parametric form", "hyperbola asymptote", "transformation rule", "curve fitting", 
    "graphing", "quadratic equations", "x-intercept", "y-intercept", "intersection point", 
    "tangent line", "conic transformation", "conic form", "rational functions",
    "length"
]
math_words.extend(coordinate_geometry_words)

# 8. 3D (Stereo) Geometry – Spatial visualization, polyhedra, spheres, tetrahedra.
stereo_geometry_terms = [
    "3D", "space", "coordinates", "point", "line", "plane", "angle", "distance", "intersection",
    "vector", "dot product", "cross product", "normal vector", "direction ratios", "spherical",
    "coordinates", "orthogonal", "parallel", "perpendicular", "plane equation", "distance from a point to a plane",
    "volume", "surface area", "tetrahedron", "cube", "sphere", "polyhedron", "dodecahedron", "icosahedron", 
    "octahedron", "prism", "pyramid", "cuboid", "cylinder", "cone", "torus", "ellipsoid", "paraboloid",
    "hyperboloid", "conical", "angle between vectors", "projection", "affine geometry", "affine transformation",
    "reflection", "rotation", "scaling", "translation", "shear", "Euler's polyhedron formula", 
    "solid angle", "dihedral angle", "diagonal", "face", "edge", "vertex", "polygon", "triangulation", 
    "adjacent faces", "opposite faces", "circumsphere", "incircle", "incenter", "circumcenter", "orthocenter", 
    "centroid", "volume formula", "surface area formula", "solid geometry", "convexity", "concavity", 
    "convex hull", "inscribed", "exscribed", "spherical geometry", "coordinate transformation", 
    "coordinates transformation", "homogeneous coordinates", "projective geometry", "stereographic projection",
    "tangent plane", "barycentric coordinates", "convex polytope", "cell complex", "simplicial complex", 
    "triangulation", "dual polyhedron", "Euler characteristic", "planar graph", "geodesic", "orbital geometry", 
    "space filling", "Klein bottle", "Möbius strip", "Lobachevskian geometry", "hyperbolic geometry", 
    "geodesic distance", "radius", "diameter", "circumference", "latitudes", "longitudes", "angular distance", 
    "Pythagorean theorem", "stereographic projection", "unit sphere", "vector cross product", 
    "tetrahedron volume", "face area", "edge length", "quadric surface", "distance formula", 
    "distance between two points", "sphere equation", "sphere center", "sphere radius", "parallel lines", 
    "slant height", "pythagorean triples", "right angle", "regular polyhedron", "semiperimeter", "diagonal plane", 
    "polar coordinates", "3D transformation", "curvature", "geodesic curvature", "planar projection", "Euler's number", 
    "spatial orientation", "octant", "coordinate plane", "cone equation", "cylinder equation", 
    "ellipsoid equation", "paraboloid equation", "spherical distance", "cylinder volume", "cylinder surface area", 
    "sphere volume", "sphere surface area", "tetrahedron surface area", "tetrahedron edges", "frustum", 
    "circumscribed sphere", "planar surface", "affine map", "linear map", "3D rotation", "affine transformation", 
    "homogeneous transformation", "projective space", "symmetry group", "equator", "parallel planes", 
    "orthogonal projection", "3D coordinates", "locus", "polar surface", "cross section", "generalized polyhedron", 
    "curved surface", "spatial reflection", "planar intersection", "line of intersection", "spatial linearity", "length"
]
math_words.extend(stereo_geometry_terms)

# 9. Trigonometry – Law of sines and cosines, identities, trigonometric equations.
trigonometry_terms = [
    "sine", "cosine", "tangent", "secant", "cosecant", "cotangent", "angle", "radian", "degree", 
    "unit circle", "trigonometric function", "identity", "law of sines", "law of cosines", "Pythagorean identity",
    "sine rule", "cosine rule", "double angle formula", "half angle formula", "sum of angles", "difference of angles", 
    "angle addition", "angle subtraction", "trigonometric equation", "trigonometric identity", "inverse sine", 
    "inverse cosine", "inverse tangent", "sine inverse", "cosine inverse", "tangent inverse", "cotangent inverse", 
    "secant inverse", "cosecant inverse", "quadrant", "reference angle", "periodicity", "amplitude", 
    "frequency", "phase shift", "harmonic", "sine wave", "cosine wave", "tangent wave", "cotangent wave", 
    "secant wave", "cosecant wave", "trigonometric ratio", "angle of elevation", "angle of depression", 
    "height", "distance", "right triangle", "acute angle", "obtuse angle", "reflex angle", "complementary angles", 
    "supplementary angles", "cotangent rule", "cosecant rule", "trig values", "trig ratios", "adjacent side", 
    "opposite side", "hypotenuse", "adjacent over hypotenuse", "opposite over hypotenuse", "opposite over adjacent", 
    "cosine law", "sine law", "angle bisector", "area of triangle", "Heron's formula", "side length", 
    "circumradius", "inradius", "trigonometric substitution", "unit circle identity", "law of tangents", 
    "secant rule", "Pythagorean theorem", "sine square plus cosine square", "trigonometric sum", "factorization", 
    "complex exponential", "Euler's formula", "circular functions", "Euler's identity", "circular angle", 
    "coordinate plane", "polar coordinates", "sine wave form", "cosine wave form", "graph of sine", "graph of cosine", 
    "phase angle", "period", "frequency", "inverse trig function", "range of sine", "range of cosine", 
    "range of tangent", "range of secant", "range of cosecant", "range of cotangent", "range of inverse trig", 
    "periodicity of sine", "periodicity of cosine", "periodicity of tangent", "periodicity of secant", 
    "periodicity of cosecant", "periodicity of cotangent", "sine function", "cosine function", "tangent function", 
    "cotangent function", "secant function", "cosecant function", "law of angles", "angle subtraction identity", 
    "addition formula", "trigonometric polynomial", "sine addition", "cosine addition", "sine difference", 
    "cosine difference", "trigonometric simplification", "trigonometric reduction", "angle multiplication formula", 
    "sine multiplication", "cosine multiplication", "double angle identity", "half angle identity", "trig function graph", 
    "trigonometric expression", "trigonometric inequality", "trigonometric ratio identity", "solution of trig equation", 
    "cotangent function", "graph of tangent", "periodicity of functions", "trig formula simplification", 
    "angle sum identity", "function period", "periodicity proof", "trigonometric series", "inverse tangent", 
    "inverse cosine", "inverse sine", "inverse secant", "inverse cosecant", "inverse cotangent", "half angle formula", 
    "law of cosines equation", "sine law formula", "tangent line", "polar angle", "angle of rotation", "right angle",
    "positive angle", "negative angle", "graph of cotangent", "circle", "circular motion", "Pythagorean identity", 
    "solution set", "unit circle equation", "circle trigonometry", "sine wave properties", "cosine wave properties"
]
math_words.extend(trigonometry_terms)

# 10. Geometric Inequalities – Triangle inequalities, Jensen’s inequality, area bounds.
geometric_inequalities_terms = [
    "triangle inequality", "Jensen's inequality", "area", "bounds", "inequality", "inequality theorem", "Heron's inequality",
    "convex function", "concave function", "mean", "convexity", "concavity", "equilateral triangle", "isosceles triangle",
    "scalene triangle", "perimeter", "side lengths", "circumradius", "inradius", "semiperimeter", "Pythagorean inequality",
    "Euler's inequality", "Cauchy-Schwarz inequality", "Minkowski inequality", "Schur inequality", "Ravi substitution", 
    "Nesbitt's inequality", "AM-GM inequality", "Cauchy inequality", "triangle", "triangle area", "geometric mean", 
    "arithmetic mean", "harmonic mean", "quadrilateral inequality", "parallelogram inequality", "circle inequality", 
    "polar inequality", "angle bisector inequality", "maximum area", "minimum area", "bounding", "bound", "upper bound", 
    "lower bound", "symmetry", "inequality chain", "Karamata's inequality", "area formula", "triangle inequality proof", 
    "equality condition", "equation", "maximum perimeter", "minimum perimeter", "radius inequality", "circumcenter", 
    "centroid", "altitude", "median", "bisector", "vertex", "side", "vertex angle", "convex polygon", "concave polygon", 
    "separation inequality", "side ratio", "perimeter bound", "height", "bounding box", "area maximization", "geometry proof",
    "quadrilateral", "parabola", "hyperbola", "ellipse", "right angle", "acute angle", "obtuse angle", "interior angle", 
    "external angle", "Pythagorean theorem", "isosceles inequality", "exterior angle inequality", "circumcircle", 
    "incircle", "descarte's rule of signs", "Hessenberg inequality", "binomial inequality", "simplex inequality", 
    "tangent", "secant", "chord", "ceva's theorem", "menelaus's theorem", "sine rule", "cosine rule", "Heron's formula", 
    "equation of lines", "equation of curves", "trigonometric inequality", "Jensen's functional inequality", 
    "lower bound theorem", "upper bound theorem", "inequality systems", "partition inequality", "convex hull", "geodesic", 
    "metric space", "convex hull area", "bounding condition", "diagonal inequality", "side length inequality", 
    "strong convexity", "triangular inequality", "convex geometry", "finite field inequality", "rearrangement inequality", 
    "fuzzy bounds", "equation constraints", "bounded set", "sublinear", "coherence", "tangent planes", "sharp inequality", 
    "Wolfram's inequality", "Bernstein's inequality", "isoperimetric inequality", "AM-HM inequality", "mean value theorem", 
    "strong convex function", "Chebyshev inequality", "unimodal function", "directed area", "Jordan curve", 
    "mathematical bound", "lower bound analysis", "quadratic inequalities", "order of magnitude", "bounding polynomial", 
    "uniform bounds", "degree of convexity", "concave hull", "finite difference inequality", "analytic inequality", 
    "subgradient", "inequality relation", "maximum area bound", "extreme values", "trivial inequality", "Lagrange's inequality"
]
math_words.extend(geometric_inequalities_terms)

## Combinatorics
# 11. Counting & Binomial Coefficients – Pigeonhole principle, stars and bars, combinatorial proofs.
counting_binomial_terms = [
    "counting", "binomial coefficient", "pigeonhole principle", "stars and bars", "combinatorial proof", "permutation",
    "combination", "factorial", "multinomial coefficient", "binomial expansion", "Pascal's triangle", "combinations with repetition",
    "derangement", "bijection", "inclusion-exclusion", "Stirling numbers", "recursion", "equivalence class", "set theory", 
    "subset", "power set", "cardinality", "partition", "bijection principle", "summation", "modular arithmetic", "integer partition",
    "combinatorial identity", "subset sum", "ordered pair", "unordered pair", "combinatorics", "order", "arrangement", "subsequence", 
    "product rule", "sum rule", "binomial identity", "principle of inclusion-exclusion", "Catalan number", "Catalan's identity",
    "Fibonacci number", "recurrence relation", "counting argument", "counting function", "multiset", "cross product", "permutation group",
    "probability", "binomial distribution", "permutation formula", "combination formula", "multinomial expansion", "counting principle", 
    "factorial identity", "combination with repetition", "counting function", "cardinal number", "maximal subset", "total number", 
    "ordered selection", "unordered selection", "subset sum problem", "fixed-point principle", "equivalence relation", "composition", 
    "generating function", "hypergeometric distribution", "stars-and-bars formula", "pigeonhole", "non-decreasing sequence", 
    "infinite sum", "geometric progression", "binomial theorem", "generalized binomial theorem", "identical objects", "repetition", 
    "consecutive elements", "ordered tuple", "unordered tuple", "proof by induction", "inductive step", "base case", "counting set", 
    "set of choices", "order of selection", "unordered choices", "combinatorics identity", "combinatorial argument", "inclusion-exclusion formula",
    "finite set", "binomial coefficient formula", "equal distribution", "equal parts", "distribution", "set cover", "counting arrangement",
    "combinatorial number", "subpermutation", "hypergeometric coefficient", "elementary counting", "permutations with repetition", 
    "combinations with repetition", "permutations of multiset", "stars-and-bars problem", "binomial expansion theorem", "counting principle",
    "counting method", "inclusion-exclusion principle", "Pólya enumeration", "graph theory", "graph isomorphism", "integer solutions", 
    "integer equations", "integer programming", "factorization", "polynomial identity", "binomial identity", "geometric sum", "factorization tree",
    "disjoint sets", "multi-step selection", "subset product", "distinct items", "sum of combinations", "recurrence formula", "combinatorial number",
    "mathematical induction", "selection problem", "combinatorial arguments", "binomial sum", "binomial coefficient properties", "dual principle",
    "k-element subset", "combinatorial number theory", "factorial decomposition", "linear recurrence", "combinatorial partition", "collection",
    "multiplicative property", "factorial expansion", "expansion formula", "expansion proof", "subset problem", "modular counting", "counting method", 
    "generalized combinations", "binomial sum formula", "counting process", "optimal selection", "elementary counting principle"
]
math_words.extend(counting_binomial_terms)

# 12. Graph Theory – Eulerian and Hamiltonian paths, trees, graph coloring, networks.
graph_theory_terms = [
    "graph", "vertex", "edge", "Eulerian path", "Hamiltonian path", "tree", "cycle", "planar graph", "connected", "degree", 
    "directed graph", "undirected graph", "bipartite graph", "subgraph", "graph coloring", "chromatic number", "adjacency", "incidence", 
    "graph traversal", "depth-first search", "breadth-first search", "Dijkstra's algorithm", "Prim's algorithm", "Kruskal's algorithm", 
    "minimum spanning tree", "strongly connected", "weakly connected", "complete graph", "isomorphic graphs", "graph isomorphism", 
    "graph density", "matching", "cut", "component", "planarity", "Euler's formula", "Hamiltonian cycle", "connected components", 
    "clique", "independent set", "vertex cover", "edge cover", "network flow", "maximum flow", "Ford-Fulkerson", "min-cut", "max-flow min-cut", 
    "k-connectivity", "graph distance", "eccentricity", "diameter", "radius", "center", "graph decomposition", "graph labeling", "reachability", 
    "distance matrix", "weighted graph", "unweighted graph", "adjacency matrix", "incidence matrix", "path", "walk", "cycle space", 
    "cut-set", "tree decomposition", "center of graph", "graph traversal algorithm", "topological sort", "bipartite matching", "graph partitioning", 
    "planarity testing", "Hamiltonian cycle problem", "Euler's theorem", "Petersen graph", "graph connectivity", "connected subgraph", "chordless cycle", 
    "graph isomorphism problem", "greedy algorithm", "minimum cut", "maximum independent set", "graph theory lemma", "coloring problem", "Four Color Theorem",
    "graph coloring algorithm", "degree sequence", "degeneracy", "spanning tree", "spanning forest", "König's theorem", "vertex-disjoint paths", "connectedness", 
    "network flow theory", "dominating set", "vertex clique", "edge-weighted graph", "tree traversal", "strongly connected components", "global minimum cut", 
    "path decomposition", "Hamiltonian path problem", "planarity condition", "degree sequence", "graph partition", "shortest path", "graph theory problem", 
    "graph theory proof", "graph theory lemma", "incidence geometry", "shortest path problem", "graph theory application", "maximal independent set", "graph diameter", 
    "matching theory", "Eulerian cycle", "generalized coloring", "maximum clique", "clique number", "connectivity function", "planar embedding", "graph algebra", 
    "signed graph", "graph analysis", "extremal graph theory", "Laplacian matrix", "graph flow", "graph traversal problem", "vertex-weighted graph", 
    "edge-weighted graph", "generalized graph", "topological graph", "non-planar graph", "graph kernel", "Hamiltonian circuit", "contraction", 
    "graph invariants", "graph dynamics", "network theory", "distance-based graph", "graph reconstruction", "radius of graph", "vertex-transitive graph", 
    "subgraph isomorphism", "graph homomorphism", "graph degree sequence", "degree sum", "algebraic graph theory", "graph decomposition theorem", "dominance relation"
]
math_words.extend(graph_theory_terms)

# 13. Invariants & Extremal Principle – Proof techniques using fixed properties and optimization.
invariants_extremal_principle_terms = [
    "invariant", "extremal principle", "optimization", "fixed point", "convexity", "concavity", "extrema", "maximum", "minimum", "local maximum", 
    "local minimum", "global maximum", "global minimum", "critical point", "Lagrange multiplier", "KKT conditions", "convex optimization", 
    "non-convex optimization", "optimization problem", "proof by contradiction", "proof by induction", "extremal problem", "optimum", "feasible solution", 
    "objective function", "constraint", "dual problem", "Karush-Kuhn-Tucker", "dual variable", "gradient", "Hessian", "gradient descent", 
    "local extremum", "global extremum", "convex set", "convex hull", "unconstrained optimization", "constrained optimization", "linearity", 
    "saddle point", "simplex method", "dual optimization", "Euler's theorem", "Rolle's theorem", "Mean Value Theorem", "Taylor series", 
    "intermediate value theorem", "function optimization", "extreme value theorem", "critical value", "local analysis", "stationary point", 
    "functional analysis", "linear programming", "integer programming", "integer optimization", "simultaneous optimization", "multi-objective optimization", 
    "quadratic optimization", "global convergence", "local convergence", "nonlinear optimization", "optimization method", "fixed point theorem", 
    "Banach fixed-point theorem", "pigeonhole principle", "sum of squares", "variational problem", "optimization theory", "dual formulation", 
    "combinatorial optimization", "graph optimization", "integer linear programming", "quadratic programming", "linear functional", "dual function", 
    "primal-dual", "optimality conditions", "bounded solution", "optimization bound", "decoupling", "separation theorem", "geometrical optimization", 
    "continuous optimization", "discrete optimization", "constraint qualification", "Kuhn-Tucker conditions", "Fejer sequence", "scaling", "singularity", 
    "subdifferential", "subgradient", "dual optimization problem", "constraint set", "relative error", "minimum value", "maximum value", 
    "optimization algorithm", "convex function", "saddle-point optimization", "maximum likelihood", "convex-concave programming", "generalized inequality", 
    "constraint matrix", "feasibility region", "monotonicity", "suboptimal solution", "optimization model", "global optimization", "linear inequality", 
    "complementary slackness", "dynamic programming", "geodesic", "quadratic form", "Lagrangian", "convex hull", "gradient ascent", "gradient descent method", 
    "sequential quadratic programming", "KKT conditions", "minimum-variance", "maximization", "duality gap", "exterior calculus", "functional derivative", 
    "gradient flow", "convergence rate", "convergence analysis", "barrier method", "cutting-plane method", "feasible point", "boundedness", "complementary slackness", 
    "dual feasible region", "coercivity", "boundary conditions", "multivariate optimization", "simplex tableau", "projection method", "primal-dual algorithm", 
    "geometry of optimization", "symmetry in optimization", "solution space", "optimization framework", "optimal solution", "bilinear form", "global solution", 
    "feasibility check", "optimum point", "extremum problem", "strong duality", "optimality theorem", "variational inequality", "maximization problem", 
    "minimization problem", "norm-based optimization", "decay rate", "convex programming", "dual optimization theory", "gradient method", "exponentially decaying", 
    "Lagrange duality", "stochastic optimization", "multi-dimensional optimization", "constraint handling", "interior-point method", "constraint function"
]
math_words.extend(invariants_extremal_principle_terms)

# 14. Recursion & Generating Functions – Counting sequences, recurrence relations, Catalan numbers.
recursion_generating_functions_terms = [
    "recursion", "recurrence relation", "generating function", "Catalan number", "sequence", "binomial coefficient", "Fibonacci sequence", "closed form", 
    "sum", "product", "partition", "counting sequence", "combinatorial", "inclusion-exclusion", "induction", "linear recurrence", "non-linear recurrence", 
    "recursive formula", "initial conditions", "step function", "polynomial", "series expansion", "recursive algorithm", "subsequence", "dividing", "factorial", 
    "infinite series", "convergence", "divisibility", "principal generating function", "ordinary generating function", "exponential generating function", 
    "recursively defined", "harmonic number", "integer partition", "asymptotic", "direct sum", "generating set", "alternating series", "divisor sum", "binomial expansion", 
    "Catalan's conjecture", "recursive tree", "summation", "hypergeometric series", "power series", "trigonometric series", "convolution", "Bernoulli number", "Eulerian number", 
    "Gamma function", "Zeta function", "recursion depth", "recursive sequence", "recursive tree", "mathematical induction", "degree", "multiplication", "counting function", 
    "monomial", "bivariate generating function", "recursion relation", "partition function", "converging", "convergent sequence", "root of equation", "binary tree", "Catalan path", 
    "combinatorial identity", "sum of powers", "generating series", "exponential generating series", "polynomial identity", "Lucas' theorem", "rational generating function", 
    "formula manipulation", "function recurrence", "specific sequence", "special functions", "Jordan's lemma", "number partition", "product expansion", "sum of divisors", 
    "interval sum", "non-homogeneous recurrence", "binomial identity", "integer solution", "combinatorial recursion", "derivative of series", "general solution", "fixed point", 
    "simplification", "Möbius inversion", "Pascal's identity", "infinite sum", "partition identity", "binomial theorem", "recursive descent", "shift operator", "multinomial theorem", 
    "power series expansion", "multiplicative function", "finite sequence", "sums and products", "recursive function", "equivalence", "upper bound", "lower bound", "telescope sum", 
    "partition number", "special sum", "asymptotic behavior", "recursive relation expansion", "recursive sequence sum", "equation solving", "nonlinear sequence", 
    "restricted partition", "recursive case", "generating function expansion", "step recursion", "generalized binomial theorem", "Laplace transform", "matrix recurrence", 
    "Euler's recurrence", "Catalan recursion", "determinants", "divisibility function", "generating tree", "combinatorial explosion", "direct sum formula", "Mathematica", 
    "orthogonal polynomials", "fractional recurrence", "recursive substitution", "continuous function", "partial sum", "summation formula", "distinct parts", "index function", 
    "recursion limit", "multiplicative recursion", "recursive back substitution", "reduction formula", "homogeneous recurrence", "maximal sequence", "non-homogeneous sum", 
    "recursively defined series", "finite difference", "shift transformation", "convergent function", "explicit formula", "sequence limit", "constant term", "catalan transformation", 
    "indexed sum", "recurrent relation", "combinatorics", "multiplicative inverse", "sum of squares", "root expansion", "symmetric function", "final answer", "root convergence", 
    "p-adic", "recursive function call", "multidimensional recurrence", "recurrent structure", "direct sum expansion", "Euler transform", "binomial recursion", "partitioning function", 
    "modular arithmetic", "combinatorics identity", "recursive terms", "multiplying sequences", "multi-step recurrence", "complexity theory", "zero solution", "extended binomial identity"
]
math_words.extend(recursion_generating_functions_terms)

# 15. Probability & Expected Value – Classic probability, linearity of expectation, Markov chains.
probability_expected_value_terms = [
    "probability", "expected value", "random variable", "distribution", "mean", "variance", "standard deviation", "probability mass function", "probability density function", 
    "Bernoulli", "binomial distribution", "uniform distribution", "normal distribution", "exponential distribution", "Poisson distribution", "Markov chain", "transition matrix", 
    "state space", "initial state", "transition probability", "absorption time", "transition rate", "conditional probability", "Bayes' theorem", "independence", "dependent events", 
    "combinatorics", "counting", "permutation", "combination", "law of total probability", "law of total expectation", "probability tree", "Markov property", "stochastic", 
    "random walk", "moment generating function", "probability theory", "independent events", "expected return", "average case", "tails", "head", "fair coin", "fair die", 
    "cumulative distribution function", "stochastic process", "joint distribution", "conditional expectation", "variance decomposition", "discrete random variable", "continuous random variable", 
    "law of large numbers", "central limit theorem", "stochastic dominance", "statistical independence", "variance-covariance", "Laplace transform", "Poisson process", "binomial coefficient", 
    "hypergeometric distribution", "normal approximation", "Chebyshev's inequality", "Chernoff bound", "Big-O notation", "expected payoff", "exponential decay", "geometric distribution", 
    "Gaussian distribution", "parametric test", "empirical probability", "probability space", "sigma-algebra", "random experiment", "sample space", "finite sample space", "infinite sample space", 
    "independent trials", "sampling without replacement", "sampling with replacement", "random sampling", "markov property", "deterministic process", "stationary distribution", "transient state", 
    "absorbing state", "expected value of a random variable", "discrete uniform", "continuous uniform", "expected time", "stationary process", "multi-step process", "randomized algorithm", 
    "success probability", "failure probability", "geometric random variable", "finite Markov chain", "steady-state distribution", "absorption probability", "mean recurrence time", "expected number of trials", 
    "conditional probability distribution", "mean first passage time", "randomized experiment", "law of iterated expectation", "Bernoulli trial", "recursive Markov chain", "convergence in probability", 
    "variance of the mean", "uniformly distributed", "biased coin", "fair dice", "probabilistic model", "regression to the mean", "probability mass", "probability density", "continuous model", 
    "expected outcome", "sample mean", "expected score", "binomial coefficient", "Markov decision process", "state transition", "path probability", "random process", "geometric mean", 
    "stochastic dominance", "entropy", "non-homogeneous Markov chain", "relative entropy", "multinomial distribution", "multivariate normal", "expected value theorem", "probabilistic bounds", 
    "momentum", "birth-death process", "polling model", "first-passage time", "probability function", "symmetry of probability", "non-homogeneous distribution", "probability distribution function", 
    "empirical distribution", "probabilistic analysis", "first moment", "second moment", "third moment", "upper bound", "lower bound", "tail probability", "multinomial coefficient", 
    "Markov property", "continuous random variable", "stochastic differential equation", "gamma distribution", "expected count", "variance of random variable", "random experiment model", 
    "Poisson distribution", "Markov decision process", "maximum likelihood estimation", "assumption of independence", "sampling distribution", "central tendency", "Poisson model", "mean of a random variable", 
    "Poisson rate", "random variable model", "empirical measure", "expected value estimator", "multivariate distribution", "point estimate", "probability estimate", "expected frequency", 
    "probability theory model", "method of moments", "bias-variance tradeoff", "linear regression model", "covariance", "sample variance", "sample mean", "probability theory", "probability density function"
]
math_words.extend(probability_expected_value_terms)

## Calculus & Advanced Topics
# 16. Limits & Continuity – Evaluating limits, squeeze theorem, L’Hôpital’s rule.
limits_continuity_terms = [
    "limit", "continuity", "squeeze theorem", "L'Hopital's rule", "derivative", "indeterminate form", "convergence", "divergence", "infinity", "epsilon", 
    "delta", "neighborhood", "left-hand limit", "right-hand limit", "two-sided limit", "bounded", "unbounded", "continuous function", "discontinuous", "jump discontinuity", 
    "removable discontinuity", "non-removable discontinuity", "point of discontinuity", "piecewise function", "asymptote", "vertical asymptote", "horizontal asymptote", 
    "oblique asymptote", "limit at infinity", "limit from the left", "limit from the right", "approaching", "infinitesimal", "singularity", "limit of a sequence", 
    "limit of a function", "approximating", "infinite limit", "finite limit", "infinite series", "convergence test", "divergent series", "continuous at a point", "uniform continuity", 
    "absolute continuity", "continuous extension", "infinitesimal behavior", "big O", "little o", "local maximum", "local minimum", "critical point", "extremum", "stationary point", 
    "gradient", "rate of change", "chain rule", "power rule", "product rule", "quotient rule", "sine rule", "cosine rule", "limit of sequence", "infinite series", 
    "convergence criteria", "Monotone Convergence Theorem", "continuity test", "differentiable", "differentiability", "differentiation", "mean value theorem", "second derivative", 
    "Taylor series", "Maclaurin series", "binomial expansion", "polynomial approximation", "growth rate", "decay rate", "L’Hopital's rule", "indeterminate form", "ratio test", 
    "limit comparison test", "root test", "squeeze theorem", "convergent sequence", "divergent sequence", "function behavior", "power series", "absolute value", "bounded function", 
    "unbounded function", "dominant term", "small-o notation", "big-O notation", "limit of derivatives", "stability", "boundedness", "limit of a product", "limit of a quotient", 
    "rule of substitution", "limit of sums", "limit of differences", "fixed point", "fixed point theorem", "uniform convergence", "uniform continuity", "pointwise convergence", 
    "limit of integrals", "limit of series", "symmetric difference", "approaching zero", "limit law", "partial limit", "limiting behavior", "limit of exponential", "limit of logarithmic", 
    "limit of trigonometric functions", "convergent integral", "infinite limits", "limit of exponential growth", "limit of logarithmic decay", "asymptotic behavior", "substitution rule", 
    "rate of convergence", "convergence of functions", "convergence of series", "convergence of sequences", "bounded variation", "function behavior at infinity", "absolute convergence", 
    "conditional convergence", "pointwise continuity", "L'Hopital's theorem", "evaluating limits", "constant sequence", "discontinuous function", "uniform convergence criteria", 
    "pointwise limit", "sequential continuity", "order of growth", "power growth", "exponential growth", "logarithmic growth", "rate of decay", "bounded sequence", "monotonicity", 
    "boundedness test", "bounded function theorem", "limit inferior", "limit superior", "infinitesimal approximation", "uniformity of limits", "continuous limit", "continuous function", 
    "local continuity", "discontinuity behavior", "singularities", "function continuity", "differentiation limits", "squeeze condition", "converging sequences", "tangent line approximation", 
    "instantaneous rate", "accelerating sequence", "limit of an infinite series", "stability of limits", "infinitesimal limits", "limit and derivative relation", "continuity condition"
]
math_words.extend(limits_continuity_terms)

# 17. Derivatives & Applications – Tangents, optimization, convexity, inequalities via derivatives.
derivatives_applications_terms = [
    "derivative", "tangent", "optimization", "convexity", "concavity", "critical point", "extrema", "local maximum", "local minimum", "global maximum", 
    "global minimum", "inflection point", "stationary point", "slope", "rate of change", "gradient", "differentiability", "mean value theorem", "chain rule", 
    "product rule", "quotient rule", "power rule", "second derivative", "first derivative", "derivative test", "monotonicity", "increasing", "decreasing", 
    "point of inflection", "concave up", "concave down", "convex function", "concave function", "convexity test", "differential", "tangent line", "normal line", 
    "optimization problem", "Lagrange multiplier", "extreme value theorem", "critical value", "function optimization", "maximum", "minimum", "optimization technique", 
    "extremum", "saddle point", "optimal solution", "cost function", "objective function", "constraint", "local optimum", "global optimum", "Euler-Lagrange equation", 
    "convex optimization", "differential equation", "differential inequality", "Taylor expansion", "Maclaurin series", "second derivative test", "continuous function", 
    "differential calculus", "implicit differentiation", "implicit function", "derivative of composite function", "higher-order derivatives", "partial derivative", "Jacobian", 
    "Hessian matrix", "directional derivative", "convex hull", "global extrema", "saddle point theorem", "Newton's method", "fixed point iteration", "monotonic sequence", 
    "mean value inequality", "derivative of inverse function", "Newton-Raphson method", "monotonic function", "critical point test", "optimization theory", "convex set", 
    "convex function properties", "decreasing function", "increasing function", "Lipschitz condition", "bounded function", "upper bound", "lower bound", "sublevel set", 
    "function behavior", "point of tangency", "tangent vector", "derivative approximation", "gradient descent", "cost optimization", "gradient ascent", "stochastic gradient", 
    "gradient of function", "derivative of logarithm", "derivative of exponential", "derivative of polynomial", "first derivative test", "second derivative test", 
    "concave function properties", "linear approximation", "convex optimization problem", "multivariable calculus", "partial derivative test", "constraint optimization", 
    "derivative of trigonometric functions", "inverse function derivative", "critical region", "second-order conditions", "derivative of rational function", 
    "optimality condition", "function minimization", "error function", "monotonic increase", "monotonic decrease", "relative maximum", "relative minimum", "differentiation rules", 
    "function increasing", "function decreasing", "boundedness", "derivative of sine", "derivative of cosine", "convexity criteria", "convexity and concavity", "inflection criterion", 
    "convex analysis", "application of derivatives", "convexity in optimization", "differential inequality", "differentiation of parametric equations", "slope of tangent", 
    "convexity of quadratic", "differential form", "derivative of hyperbolic functions", "extremum theorem", "local analysis", "global analysis", "derivative of square root", 
    "concave shape", "convex shape", "constraint function", "derivatives of power functions", "function analysis", "differentiability of function", "Lagrange multiplier method", 
    "tangent approximation", "optimization with constraints", "first order condition", "second order condition", "derivative of arctangent", "derivative of arccosine", 
    "second order derivative", "function optimization problem", "partial derivative test", "derivative of inverse trigonometric", "critical value calculation", "optimization strategy", 
    "convex optimization theory", "differential constraint", "extremum points", "convexity in optimization", "differentiation formula", "concavity criteria", "second derivative approximation"
]
math_words.extend(derivatives_applications_terms)

# 18. Integral Calculus – Summation approximations, definite and indefinite integrals.
integral_calculus_terms = [
    "integral", "summation", "approximation", "definite integral", "indefinite integral", "integration", "antiderivative", "integral bounds", "Riemann sum", 
    "Lebesgue integral", "integrand", "integration by parts", "substitution", "u-substitution", "integration by partial fractions", "integral formula", "improper integral", 
    "convergence", "divergence", "area under curve", "fundamental theorem of calculus", "mean value theorem", "integral approximation", "integral convergence", "differential", 
    "continuous function", "integral identity", "integral method", "integral table", "integral limits", "integral solution", "curve integral", "line integral", 
    "surface integral", "volume integral", "double integral", "triple integral", "multiple integral", "polar coordinates", "cartesian coordinates", "spherical coordinates", 
    "change of variables", "integral evaluation", "definite limit", "indefinite limit", "integration constant", "constant of integration", "antiderivative of function", 
    "inverse function theorem", "quadrature", "midpoint rule", "trapezoidal rule", "Simpson's rule", "integration technique", "reduction formula", "derivative", "differentiation", 
    "function limits", "Riemann integral", "Cauchy principal value", "Euler's method", "numerical integration", "adaptive quadrature", "Gaussian quadrature", "Maclaurin series", 
    "Taylor expansion", "approximate area", "power series", "series expansion", "beta function", "gamma function", "integral transform", "Laplace transform", "Fourier transform", 
    "infinite series", "convergence test", "divergent series", "power rule", "integral convergence test", "partial fractions", "linear operator", "L'Hopital's rule", "symmetric integrals", 
    "sum of series", "sum of terms", "sum of integrals", "rectangular approximation", "approximation of integrals", "integral bounds", "continuity", "smooth function", 
    "partial integration", "integral representation", "inverse Laplace transform", "Laplace transform properties", "Fourier series", "Fourier integral", "Fourier transform", 
    "integral equation", "solution by integration", "reducing integral", "integration by parts formula", "definite area", "Cauchy integral theorem", "contour integral", 
    "generalized integral", "non-elementary integrals", "Riemann-Stieltjes integral", "Green's theorem", "Stokes' theorem", "divergence theorem", "calculus of variations", 
    "integral of exponential", "integral of logarithm", "integral of trigonometric", "integral of inverse trigonometric", "integral of hyperbolic", "integral of rational functions", 
    "integration of polynomial", "rational function integration", "Laplace's method", "asymptotic approximation", "series expansion", "differential equations", "integral equation solution", 
    "general integral", "convergence of integrals", "analytic function", "Hermite polynomial", "Legendre polynomial", "integral limits evaluation", "arbitrary constant", 
    "variable substitution", "multivariable integral", "double integration", "triple integration", "integral constraints", "integral properties", "integral bounds comparison", 
    "finite sum", "infinite sum", "integral sum", "factorization", "reduction of integral", "integration of rational functions", "singular integral", "Stokes' theorem proof", 
    "Taylor series expansion", "polar integral", "polar form", "hyperbolic functions", "boundary conditions", "solving integrals", "indeterminate forms", "quadrature formulas", 
    "integral region", "principal value", "multiple summation", "integral function", "residue theorem", "integration method", "multiple integral evaluation", "beta function integral", 
    "integral in polar coordinates", "error bounds in integration", "intermediate value theorem", "integral method of solution", "continuous differentiable", "analytic integral", 
    "complex integration", "integration of power functions", "bounded integral", "improper integral convergence", "integration of rational polynomials", "contour integration", 
    "non-integrable functions", "boundary value problem", "integral test", "integral approximation methods"
]
math_words.extend(integral_calculus_terms)

# 19. Functional Equations – Identifying functions satisfying given conditions, Cauchy’s FE, Jensen’s FE.
functional_equations_terms = [
    "functional equation", "function", "identity", "Cauchy's functional equation", "Jensen's functional equation", "additive function", "multiplicative function", 
    "homogeneous function", "functional inequality", "linear function", "quadratic function", "continuous function", "solution", "domain", "range", "mapping", "constant function", 
    "piecewise function", "exponential function", "polynomial function", "periodic function", "trigonometric function", "linear solution", "non-linear solution", 
    "Cauchy condition", "Jensen condition", "additivity", "multiplicativity", "functional form", "equation property", "constraint", "equation solver", "independence", 
    "boundedness", "discrete solution", "continuous solution", "real-valued function", "vector space", "invertible function", "functional identity", "iterative solution", 
    "set of functions", "fixed point", "monotonicity", "piecewise definition", "convex function", "concave function", "linear operator", "functional derivation", "solution uniqueness", 
    "general solution", "consequence", "transform", "expansion", "system of equations", "equation symmetry", "composition", "bijective function", "inverse function", 
    "transformation", "substitution", "generalization", "equation structure", "convergence", "congruence", "mapping property", "function property", "quadratic form", 
    "equation invariant", "function consistency", "proof by induction", "solution set", "bounded function", "recursive solution", "functional formulation", "symmetry of solution", 
    "relation", "system of functional equations", "method of solving", "equation consistency", "function rule", "complexity", "differential equation", "functional derivation", 
    "sum of functions", "product of functions", "factorization", "additive property", "multiplicative property", "pair of functions", "equation solution", "inverse relation", 
    "equation transformation", "function behavior", "sequence of functions", "bounds of function", "approximation", "analytic solution", "continuous mapping", "piecewise solution", 
    "interval", "mapping rule", "analysis", "non-linear equation", "function transformation", "solution construction", "derivation", "solution proof", "value of function", 
    "zero of function", "asymptotic behavior", "finite solution", "non-zero solution", "equation modeling", "solution existence", "monotone function", "approximation method", 
    "absolute value", "solution uniqueness proof", "bounded domain", "infinite solution", "function constant", "scalar function", "additive property proof", "continuous transformation", 
    "closure property", "function boundary", "non-decreasing function", "inverse of function", "solution by substitution", "function identity proof", "power function", "bounded equation", 
    "substitution method", "solution structure", "non-constant function", "non-linear behavior", "solution consistency", "recursive method", "Cauchy functional equation", "Cauchy’s condition", 
    "equation parameters", "transformation rule", "value theorem", "sequence solution", "discrete solution", "equation mapping", "parametric solution", "quadratic solution", 
    "continuity property", "periodicity", "function symmetry", "solution interval", "solution existence proof", "recursive equation", "equation uniqueness", "Cauchy’s proof", 
    "additive solution", "invariant function", "congruence condition", "bound of equation", "finite solution set", "recursive structure", "symmetric function", "transformation proof"
]
math_words.extend(math_words)

# 20. Game Theory & Strategy – Winning strategies, Nim games, Sprague-Grundy theorem.
game_theory_terms = [
    "game theory", "strategy", "winning strategy", "Nim game", "Sprague-Grundy theorem", "combinatorial game", "game tree", "minimax", "zero-sum game", "payoff matrix", 
    "game value", "Nash equilibrium", "subgame perfect equilibrium", "dominant strategy", "mixed strategy", "pure strategy", "two-player game", "multi-player game", 
    "optimal strategy", "best response", "prisoner's dilemma", "chess", "tic-tac-toe", "checkers", "win condition", "losing position", "Grundy number", "Nimber", 
    "impartial game", "normal play", "misère play", "position analysis", "Zermelo's theorem", "minimax theorem", "Pareto optimality", "pure strategy Nash equilibrium", 
    "mixed strategy Nash equilibrium", "extensive form", "normal form", "repeated game", "stochastic game", "best strategy", "game outcome", "decision tree", 
    "payoff function", "regret", "dominance", "sequential game", "simultaneous game", "subgame", "cooperative game", "non-cooperative game", "auction theory", 
    "auction game", "two-player zero-sum", "generalized game", "game analysis", "cooperative strategy", "non-cooperative strategy", "value iteration", "backward induction", 
    "forward induction", "symmetric game", "asymmetric game", "public goods game", "utility function", "voting game", "correlation", "commitment", "information asymmetry", 
    "bluffing", "incomplete information", "bounded rationality", "heuristics", "punishment strategy", "reward strategy", "repeated dilemma", "collaborative game", 
    "competitive game", "differential game", "auction theory", "Stag hunt", "coordination game", "sequential moves", "simultaneous moves", "strategic interaction", 
    "game dynamics", "repeated play", "move sequence", "opponent's move", "theory of moves", "game solution", "strategy profile", "strategy space", "equilibrium strategy", 
    "winning move", "losing move", "minimax strategy", "maximin strategy", "randomized strategy", "dominance relation", "payoff matrix analysis", "expected payoff", 
    "optimal choice", "backward induction", "decision making", "move set", "game partition", "Nim-sum", "Sprague-Grundy number", "equilibrium point", "minimax algorithm", 
    "multi-stage game", "game strategy", "game outcome analysis", "suboptimal play", "perfect information", "imperfect information", "information set", "discount factor", 
    "strategy iteration", "outcome matrix", "self-interested agent", "equilibrium point", "local optimum", "global optimum", "deviation strategy", "potential game", 
    "social choice theory", "game behavior", "game matrix", "payoff distribution", "solving game", "solution concept", "strategy space exploration", "strategic thinking", 
    "turn-based game", "time-dependent strategy", "block game", "exclusionary game", "dominance strategy", "coalition game", "survival game", "game payoff", 
    "supergame", "convergence theorem", "simultaneous strategy", "mixing strategy", "game classification", "game model", "sequential decision", "sequential equilibrium", 
    "game prediction", "solution path", "optimal play", "game outcome prediction", "path analysis", "neutrally optimal strategy", "mixed Nash equilibrium", "game equilibrium", 
    "strategic decision", "game practice", "game optimization", "probabilistic strategy", "exploration vs exploitation", "risk dominance", "game setup", "game space"
]
math_words.extend(game_theory_terms)

math_words = list(set(math_words))