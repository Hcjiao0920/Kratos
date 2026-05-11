// NEGATIVE-TEST FIXTURE. This file MUST FAIL TO COMPILE.
// Do not edit it to make it valid; if you need a working test, create
// a new one. This file exists only so that
// `verify_build_failure_propagates.bat` can confirm that
// `build_and_run.bat` correctly reports a non-zero exit code when the
// compiler fails (Codex review finding: a stale executable left on
// disk could otherwise mask a failed compile).

#include "../sand_hypoplastic_kernel.h"

int main() {
    // Reference an identifier that does not exist in any header. The
    // compiler must reject this with an error, not a warning.
    deliberately_undeclared_identifier_for_negative_build_test();
    return 0;
}
