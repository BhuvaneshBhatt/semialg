from semialg import errors, exceptions


def test_exception_module_reexports_canonical_classes():
    assert exceptions.SemialgError is errors.SemialgError
    assert exceptions.BackendFailure is errors.BackendFailure
    assert exceptions.UnsupportedFragmentError is errors.UnsupportedFragmentError
    assert exceptions.CertificationFailure is errors.CertificationFailure
    assert exceptions.ResourceLimitError is errors.ResourceLimitError
