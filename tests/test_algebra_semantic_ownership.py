"""Architecture regressions for semantic ownership in exact algebra modules."""

from semialg.algebraic import border_basis, gtz_primary, modular, roots
from semialg.algebraic._border_basis_types import BorderBasisDiagnostics
from semialg.algebraic._gtz_primary_types import GTZPrimaryDecompositionResult
from semialg.algebraic._modular_types import ModularGroebnerResult
from semialg.algebraic._root_certificate_types import PolynomialRootIntervalCertificate


def test_public_algorithm_modules_reexport_owned_result_models():
    assert border_basis.BorderBasisDiagnostics is BorderBasisDiagnostics
    assert gtz_primary.GTZPrimaryDecompositionResult is GTZPrimaryDecompositionResult
    assert modular.ModularGroebnerResult is ModularGroebnerResult
    assert roots.PolynomialRootIntervalCertificate is PolynomialRootIntervalCertificate


def test_hot_algorithm_entry_points_remain_directly_owned():
    assert modular.modular_groebner_basis_qq.__module__ == "semialg.algebraic.modular"
    assert gtz_primary.gtz_primary_decomposition.__module__ == "semialg.algebraic.gtz_primary"
    assert roots.isolate_real_roots.__module__ == "semialg.algebraic.roots"
    assert border_basis.compute_border_basis.__module__ == "semialg.algebraic.border_basis"


def test_certificate_replay_and_backends_have_distinct_semantic_owners():
    import semialg.algebraic_function_fields as function_fields

    assert (
        function_fields.certified_factor_univariate.__module__
        == "semialg._function_field_factorization"
    )
    assert (
        function_fields.compress_primitive_element.__module__
        == "semialg._function_field_compression"
    )
    assert (
        modular.verify_modular_groebner_certificate.__module__
        == "semialg.algebraic._modular_certification"
    )
    assert modular.modular_resultant_qq.__module__ == "semialg.algebraic._modular_resultants"
    assert (
        gtz_primary.verify_gtz_primary_decomposition_certificate.__module__
        == "semialg.algebraic._gtz_primary_certification"
    )
    assert (
        roots.certify_polynomial_root_interval.__module__ == "semialg.algebraic._root_certification"
    )
