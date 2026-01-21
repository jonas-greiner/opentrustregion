# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import unittest

# load standard tests
from pyopentrustregion.tests import (
    OpenTrustRegionTests,
    CInterfaceTests,
    PyInterfaceTests,
    SystemTests,
)

# try to load extension module tests if available
try:
    from pyopentrustregion.extensions.quasi_newton.tests import (
        QNCInterfaceTests,
        QNPyInterfaceTests,
    )
except AttributeError:
    pass

try:
    from pyopentrustregion.extensions.arh.tests import (
        ARHCInterfaceTests,
        ARHPyInterfaceTests,
    )
except AttributeError:
    pass

try:
    from pyopentrustregion.extensions.s_gek.tests import (
        SGEKCInterfaceTests,
        SGEKPyInterfaceTests,
    )
except AttributeError:
    pass


if __name__ == "__main__":
    unittest.main(verbosity=0)
