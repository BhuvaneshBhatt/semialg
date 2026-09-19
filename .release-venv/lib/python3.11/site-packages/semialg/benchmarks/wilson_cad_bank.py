"""Wilson CAD Example Bank v4 converted to semialg benchmark records."""

from __future__ import annotations

from .cases import BenchmarkCase, BenchmarkDomain, Difficulty, Provenance

SOURCE = Provenance(
    citation="Bradford, Davenport & Wilson, A Repository for CAD Examples; CAD Example Bank v4",
    url="https://doi.org/10.15125/BATH-00069",
    identifier="CAD Example Bank v4",
    license="CC BY-SA 4.0",
)

# Entries preserve the source procedure/index identifier and its prescribed variable order.
# Polynomial strings use Python exponentiation but otherwise retain the mathematical input.
_RAW = [
    ("CMXYExamples:1", "Parametric Parabola", ("a*x**2+b*x+c",), ("x", "c", "b", "a")),
    ("CMXYExamples:2", "Whitney Umbrella", ("x-u*v", "y-v", "z-u**2"), ("v", "u", "z", "y", "x")),
    ("CMXYExamples:3", "Quartic", ("x**4+p*x**2+q*x+r",), ("x", "p", "q", "r")),
    (
        "CMXYExamples:4",
        "Sphere and Catastrophe",
        ("z**2+y**2+x**2-1", "z**3+x*z+y"),
        ("x", "y", "z"),
    ),
    ("CMXYExamples:5", "Arnon-84", ("y**4-2*y**3+y**2-3*x**2*y+2*x**4",), ("y", "x")),
    (
        "CMXYExamples:6",
        "Arnon-84-2",
        ("144*y**2+96*x**2*y+9*x**4+105*x**2+70*x-98", "x*y**2+6*x*y+x**3+9*x"),
        ("y", "x"),
    ),
    (
        "CMXYExamples:7",
        "A real implicitization problem",
        ("x-u*v", "y-u*v**2", "z-u**2"),
        ("v", "u", "z", "y", "x"),
    ),
    (
        "CMXYExamples:8",
        "Ball and Circular Cylinder",
        ("x**2+y**2+z**2-1", "x**2+(y+z-2)**2-1"),
        ("z", "y", "x"),
    ),
    (
        "CMXYExamples:9",
        "Termination of term rewrite system",
        ("x-r", "y-r", "x**2*(1+2*y)**2-y**2*(1+2*x**2)"),
        ("y", "x", "r"),
    ),
    (
        "CMXYExamples:10",
        "Collins and Johnson",
        (
            "3*a**2*r+3*b**2-2*a*r-a**2-b**2",
            "3*a**2*r+3*b**2*r-4*a*r+r-2*a**2-2*b**2+2*a",
            "a-1/2",
            "b",
            "r",
            "r-1",
        ),
        ("r", "a", "b"),
    ),
    (
        "CMXYExamples:11",
        "Range of lower bounds",
        ("a", "a*z**2+b*z+c", "a*x**2+b*x+c-y"),
        ("z", "c", "b", "a", "x", "y"),
    ),
    (
        "CMXYExamples:12",
        "X-axis ellipse problem",
        ("b**2*(x-c)**2+a**2*y**2-a**2*b**2", "x**2+y**2-1"),
        ("y", "x", "b", "c", "a"),
    ),
    (
        "CMXYExamples:13",
        "Davenport and Heintz",
        ("a-d", "b-c", "a-c", "b-1", "a**2-b"),
        ("a", "b", "c", "d"),
    ),
    (
        "CMXYExamples:14",
        "Hong-90",
        ("r+s+t", "r*s+s*t+t*r-a", "r*s*t-b"),
        ("b", "a", "t", "s", "r"),
    ),
    (
        "CMXYExamples:15",
        "Solotareff-3",
        (
            "r",
            "r-1",
            "u+1",
            "u-v",
            "v-1",
            "3*u**2+2*r*u-a",
            "3*v**2+2*r*v-a",
            "u**3+r*u**2-a*u+a-r-1",
            "v**3+r*v**2-a*v-2*b-a+r+1",
        ),
        ("v", "u", "b", "r", "a"),
    ),
    (
        "CMXYExamples:16",
        "Collision Problem",
        (
            "17*t/16-6",
            "17*t/16-10",
            "x-17*t/16+1",
            "x-17*t/16-1",
            "y-17*t/16+9",
            "y-17*t/16+7",
            "(x-t)**2+y**2-1",
        ),
        ("t", "x", "y"),
    ),
    (
        "CMXYExamples:17",
        "McCallum Trivariate Random Polynomial",
        ("(y-1)*z**4+x*z**3+x*(1-y)*z**2+(y-x-1)*z+y",),
        ("z", "y", "x"),
    ),
    (
        "CMXYExamples:18",
        "Ellipse Problem",
        ("b**2*(x-c)**2+a**2*(y-d)**2-a**2*b**2", "a", "b", "x**2+y**2-1"),
        ("y", "x", "d", "c", "b", "a"),
    ),
]

# Remaining v4 families. The source bank itself does not attach canonical cell counts to
# every entry; reference counts are therefore stored only when supplied by a cited experiment.
_RAW += [
    (
        "BranchExamples:1",
        "Kahan A",
        (
            "4*y*(2*y**2*x+2*x**3+5*y**2+21*x**2+72*x+81)",
            "-225*x**2-324*x+63*y**2-4*x**4-52*x**3+12*y**2*x+4*y**4",
            "2*y",
            "2*x+9",
            "8*y",
            "8*x**2+56*x+8*y**2+96",
            "y",
            "x**2+7*x+y**2+12",
            "4*y*(2*y**2*x+2*x**3+5*y**2+21*x**2+72*x+81)",
            "-4*x**4-52*x**3-252*x**2+12*y**2*x-540*x+36*y**2+4*y**4-432",
            "4*x**4+52*x**3+225*x**2-12*y**2*x+324*x-63*y**2-4*y**4",
            "2*y",
            "-2*x-6",
            "2*x",
            "8*y",
            "-8*x**2-56*x-8*y**2-96",
            "2*x**2+8*x+2*y**2",
        ),
        ("y", "x"),
    ),
    ("BranchExamples:2", "Kahan B", (), ("x", "y")),
    (
        "BranchExamples:3",
        "ArcSin A",
        (
            "-2*x*y",
            "y**2-x**2+1",
            "-8*x*y*(2*x**2-2*y**2-1)",
            "1-4*x**2+4*y**2+4*x**4-24*x**2*y**2+4*y**4",
            "y",
            "1-x",
            "y",
            "x+1",
        ),
        ("y", "x"),
    ),
    ("BranchExamples:4", "ArcSin B", (), ("x", "y")),
    (
        "MotionExamples:1",
        "Piano Mover's Problem (Davenport)",
        (
            "(x-a)**2+(y-b)**2-9",
            "y*b",
            "x*(y-b)**2+y*(a-x)*(y-b)",
            "(y-1)*(b-1)",
            "(x+1)*(y-b)**2+(y-1)*(a-x)*(y-b)",
            "x*a",
            "y*(x-a)**2+x*(b-y)*(x-a)",
            "(x+1)*(a+1)",
            "(y-1)*(x-a)**2+(x+1)*(b-y)*(x-a)",
        ),
        ("x", "y", "a", "b"),
    ),
]
# paired entries reuse the same polynomial list with the alternate prescribed order
_RAW[-4] = (_RAW[-4][0], _RAW[-4][1], _RAW[-5][2], _RAW[-4][3])
_RAW[-2] = (_RAW[-2][0], _RAW[-2][1], _RAW[-3][2], _RAW[-2][3])

_BH_POLYS = [
    ("x**2-(1/2)*(y**2)-(1/2)*z**2", "x*z+z*y-2*x", "z**2-y"),
    ("4*x**2+x*y**2-z+(1/4)", "2*x+y**2*z+(1/2)", "x**2*z-(1/2)*x-y**2"),
    ("x**2+y**2-1", "b**2*(x-c)**2+a**2*y**2-a**2*b**2", "a", "a-1", "b", "b-1", "c", "c-1"),
    (
        "3*x**2-2*x-a",
        "x**3-x**2-a*x-2*b+a-2",
        "3*y**2-2*y-a",
        "y**3-y**2-a*y-a+2",
        "4*a-1",
        "4*a-7",
        "4*b+3",
        "4*b-3",
        "x+1",
        "x",
        "y",
        "y-1",
    ),
    ("(1/4)*(x-t)**2+(y-10)**2-1", "(1/4)*(x-a*t)**2+(y-a*t)**2-1", "t", "a"),
]
for i, (name, poly, orders) in enumerate(
    [
        ("Intersection", _BH_POLYS[0], [("z", "y", "x"), ("z", "x", "y")]),
        ("Random", _BH_POLYS[1], [("z", "y", "x"), ("x", "y", "z")]),
        ("Ellipse", _BH_POLYS[2], [("y", "x", "c", "b", "a"), ("x", "y", "c", "b", "a")]),
        ("Solotareff", _BH_POLYS[3], [("y", "x", "b", "a"), ("y", "x", "a", "b")]),
        ("Collision", _BH_POLYS[4], [("y", "x", "t", "a"), ("t", "x", "y", "a")]),
    ]
):
    base = 2 * i + 1
    _RAW.extend(
        [
            (f"BHExamples:{base}", f"{name} A", poly, orders[0]),
            (f"BHExamples:{base + 1}", f"{name} B", poly, orders[1]),
        ]
    )

_OTHER = [
    (
        "Off-Center Ellipse",
        ("a", "16*a**2*y**2-8*a**2*y+4*x**2-4*x-3*a**2+1", "y**2+x**2-1"),
        ("y", "x", "a"),
    ),
    ("Concentric Circles", ("x**2+y**2-9", "x**2+y**2-1"), ("y", "x")),
    ("Non-Concentric Circles", ("x**2+y**2-9", "x**2+(y-1)**2-1"), ("y", "x")),
    (
        "Edges Square Product",
        (
            "x-a*b+c",
            "y-a*c-x-2",
            "a",
            "2-a",
            "b-2",
            "4-b",
            "c+1",
            "1-c",
            "x+1",
            "9-x",
            "y+6",
            "6-y",
        ),
        ("c", "b", "a", "y", "x"),
    ),
    (
        "Simplified Edges Square Product",
        (
            "x-1",
            "2-a",
            "1+x",
            "9-x",
            "y+6",
            "6-y",
            "-(a**2+1)*(-y-a*x+2*a**2+2)",
            "(a**2+1)**2-(-x+a*y)*(a**2+1)",
            "a**2+1",
            "(a**2+1)*(a**2+1-x+a*y)",
            "4*(a**2+1)**2-(a**2+1)*(y+a*x)",
        ),
        ("a", "y", "x"),
    ),
    (
        "Putnum Example",
        ("a**2+b**2-1", "(c-10)**2+d**2-9", "2*x-(a+c)", "2+y-(b+d)"),
        ("d", "c", "b", "a", "y", "x"),
    ),
    (
        "Simplified Putnum Example",
        ("a**2+4*y**2-4*y*d+d**2-1", "4*x**2-4*x*a-40*x+a**2+20*a+91+d**2"),
        ("d", "a", "y", "x"),
    ),
    (
        "YangXia",
        (
            "a**2*h**2-4*s*(s-a)*(s-b)*(s-c)",
            "2*R*h-b*c",
            "2*s-a-b-c",
            "b",
            "c",
            "R",
            "h",
            "a+b-c",
            "b+c-a",
            "c+a-b",
        ),
        ("c", "b", "s", "R", "h", "a"),
    ),
    (
        "Simplified YangXia",
        (
            "-(1/2)*b",
            "R",
            "b",
            "h",
            "(1/16)*a**2*h**2*b**4-(1/32)*a**2*b**6-(1/8)*a**2*R**2*h**2*b**2-(1/8)*R**2*h**2*b**4+(1/64)*b**8+(1/64)*a**4*b**4+(1/4)*R**4*h**4",
            "-(1/4)*(-a*b-b**2+2*R*h)*b",
            "(1/2)*R*h*b",
            "(1/4)*(2*R*h+a*b-b**2)*b",
            "(1/4)*(b**2+2*R*h-a*b)*b",
        ),
        ("b", "a", "h", "R"),
    ),
    (
        "SEIT Model",
        (
            "d-d*s-b*J*s",
            "v*F-(d+t)*J",
            "b*J+c*J*T-(d+v+r)*F+(1-q)*t*J",
            "-d*T+r*F+q*t*J-b**2*T*J",
            "F",
            "J",
            "T",
            "s",
            "b",
            "d",
            "v",
            "r",
            "t",
            "q",
            "b-c",
        ),
        ("T", "J", "F", "s", "v", "t", "r", "q", "d", "c", "b"),
    ),
    (
        "Simplified SEIT Model",
        (
            "d",
            "r",
            "t",
            "q",
            "b-c",
            "v",
            "J",
            "b",
            "c",
            "d+J*b",
            "-v",
            "(d+t)*J*v",
            "v*c",
            "d*(d+J*b)",
        ),
        ("J", "v", "t", "r", "q", "d", "c", "b"),
    ),
    ("Cyclic-3", ("a+b+c", "a*b+b*c+c*a", "a*b*c-1"), ("a", "b", "c")),
    (
        "Cyclic-4",
        ("a+b+c+d", "a*b+b*c+c*d+d*a", "a*b*c+b*c*d+c*d*a+d*a*b", "a*b*c*d-1"),
        ("a", "b", "c", "d"),
    ),
    (
        "Joukowsky Transformation",
        (
            "a*(c**2+d**2)*(a**2+b**2+1)-c*(a**2+b**2)*(c**2+d**2+1)",
            "b*(c**2+d**2)*(a**2+b**2-1)-d*(a**2+b**2)*(c**2+d**2-1)",
            "b*d",
            "c**2+d**2-1",
            "a-c",
            "b-d",
        ),
        ("d", "c", "b", "a"),
    ),
]
for i, (name, polys, order) in enumerate(_OTHER, 1):
    _RAW.append((f"OtherExamples:{i}", name, polys, order))
_J = _OTHER[-1][1]
_RAW += [
    ("JoukowskyTransformation:1", "Joukowsky Transformation", _J, ("d", "c", "b", "a")),
    ("JoukowskyTransformation:2", "Separate Clauses", _J[:-1], ("d", "c", "b", "a")),
    (
        "JoukowskyTransformation:3",
        "Upper Half Plane",
        (_J[0], _J[1], "b", "d", "a-c", "b-d"),
        ("d", "c", "b", "a"),
    ),
]

# The 18 TTICAD worked examples are also part of Example Bank v4. Their formula
# structure is represented in tticad_dataset.py; include their source identities here.
_TTI_NAMES = (
    "Intersection dagger A",
    "Intersection dagger B",
    "Random dagger A",
    "Random dagger B",
    "Ellipse dagger A",
    "Ellipse dagger B",
    "Solotareff dagger A",
    "Solotareff dagger B",
    "Collision dagger A",
    "Collision dagger B",
    "Kahan A",
    "Kahan B",
    "ArcSin A",
    "ArcSin B",
    "2D Example A",
    "2D Example B",
    "3D Example A",
    "3D Example B",
)


def _base_cases() -> tuple[BenchmarkCase, ...]:
    cases = []
    for source_id, name, polys, order in _RAW:
        cases.append(
            BenchmarkCase(
                name=(
                    "joukowsky_transformation_1"
                    if source_id == "JoukowskyTransformation:1"
                    else "davenport_heintz"
                    if name == "Davenport and Heintz"
                    else (
                        "wilson_"
                        + name.lower().replace(" ", "_").replace("'", "").replace("-", "_")
                        if name in {"Cyclic-3", "Joukowsky Transformation"}
                        else name.lower().replace(" ", "_").replace("'", "").replace("-", "_")
                    )
                ),
                domain=BenchmarkDomain.CAD,
                description=name,
                variables=order,
                polynomials=polys,
                tags=("literature", "wilson-cad-bank"),
                provenance=Provenance(SOURCE.citation, SOURCE.url, source_id, SOURCE.license),
                difficulty=Difficulty(len(order), polynomial_count=len(polys)),
                expected={"source_identifier": source_id, "variable_order": order},
            )
        )
    # TTICAD source entries are supplied by the richer TTICAD records, avoiding duplicate formulas.
    return tuple(cases)


def wilson_bank_source_identifiers() -> tuple[str, ...]:
    return tuple(x[0] for x in _RAW) + tuple(f"TTICADExamples:{i}" for i in range(1, 19))


# Complete the v4 bank with its 18 TTICAD worked-example polynomial sets.
_TTICAD_RAW = [
    ("Intersection dagger A", ("x*z+z*y-2*x", "z**2-y", "2*x**2-y**2-z**2"), ("z", "y", "x")),
    ("Intersection dagger B", ("x*z+z*y-2*x", "z**2-y", "2*x**2-y**2-z**2"), ("z", "x", "y")),
    (
        "Random dagger A",
        ("16*x**2+4*x*y**2-4*z+1", "4*x+2*y**2*z+1", "2*x**2*z-x-2*y**2"),
        ("z", "y", "x"),
    ),
    (
        "Random dagger B",
        ("16*x**2+4*x*y**2-4*z+1", "4*x+2*y**2*z+1", "2*x**2*z-x-2*y**2"),
        ("z", "x", "y"),
    ),
    ("Ellipse dagger A", _BH_POLYS[2], ("y", "x", "c", "b", "a")),
    ("Ellipse dagger B", _BH_POLYS[2], ("x", "y", "c", "b", "a")),
    ("Solotareff dagger A", _BH_POLYS[3], ("y", "x", "b", "a")),
    ("Solotareff dagger B", _BH_POLYS[3], ("y", "x", "a", "b")),
    (
        "Collision dagger A",
        ("x**2-2*x*t+t**2+4*y**2-80*y+396", "t", "x**2-2*x*a*t+5*a**2*t**2+4*y**2-8*y*a*t-4", "a"),
        ("y", "x", "t", "a"),
    ),
    (
        "Collision dagger B",
        ("x**2-2*x*t+t**2+4*y**2-80*y+396", "t", "x**2-2*x*a*t+5*a**2*t**2+4*y**2-8*y*a*t-4", "a"),
        ("t", "x", "y", "a"),
    ),
    ("Kahan A", _RAW[18][2], ("y", "x")),
    ("Kahan B", _RAW[18][2], ("x", "y")),
    (
        "ArcSin A",
        (
            "-2*x*y",
            "y**2-x**2+1",
            "16*y**3*x-16*y*x**3+8*y*x",
            "1-4*x**2+4*y**2+4*x**4-24*x**2*y**2+4*y**4",
            "y",
            "1-x",
            "y",
            "x+1",
        ),
        ("y", "x"),
    ),
    (
        "ArcSin B",
        (
            "-2*x*y",
            "y**2-x**2+1",
            "16*y**3*x-16*y*x**3+8*y*x",
            "1-4*x**2+4*y**2+4*x**4-24*x**2*y**2+4*y**4",
            "y",
            "1-x",
            "y",
            "x+1",
        ),
        ("x", "y"),
    ),
    (
        "2D Example A",
        ("x**2+y**2-1", "4*x*y-1", "x**2-8*x+16+y**2-2*y", "4*x*y-4*x-17"),
        ("y", "x"),
    ),
    (
        "2D Example B",
        ("x**2+y**2-1", "4*x*y-1", "x**2-8*x+16+y**2-2*y", "4*x*y-4*x-17"),
        ("x", "y"),
    ),
    (
        "3D Example A",
        (
            "x**2+y**2+z**2-1",
            "4*x*y*z-1",
            "x**2-8*x+y**2-2*y+z**2-4*z+20",
            "4*y*x*z-8*y*x-4*x*z+8*x-16*y*z+32*y+16*z-33",
        ),
        ("z", "y", "x"),
    ),
    (
        "3D Example B",
        (
            "x**2+y**2+z**2-1",
            "4*x*y*z-1",
            "x**2-8*x+y**2-2*y+z**2-4*z+20",
            "4*y*x*z-8*y*x-4*x*z+8*x-16*y*z+32*y+16*z-33",
        ),
        ("z", "x", "y"),
    ),
]


def wilson_cad_cases() -> tuple[BenchmarkCase, ...]:
    cases = list(_base_cases())
    for i, (name, polys, order) in enumerate(_TTICAD_RAW, 1):
        source_id = f"TTICADExamples:{i}"
        cases.append(
            BenchmarkCase(
                name=f"wilson_tticad_{i:02d}",
                domain=BenchmarkDomain.CAD,
                description=name,
                variables=order,
                polynomials=polys,
                tags=("literature", "wilson-cad-bank", "tticad-worked-example"),
                provenance=Provenance(SOURCE.citation, SOURCE.url, source_id, SOURCE.license),
                difficulty=Difficulty(len(order), polynomial_count=len(polys)),
                expected={"source_identifier": source_id, "variable_order": order},
            )
        )
    return tuple(cases)
