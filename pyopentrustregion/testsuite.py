# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import sys
import unittest
from pathlib import Path

# when run as a script, import the package this file belongs to rather than another
# one that is installed
if not __package__:
    sys.path.insert(0, str(Path(__file__).parent.parent.absolute()))

# load standard tests
from pyopentrustregion.tests import (
    OpenTrustRegionUnitTests,
    CInterfaceUnitTests,
    PyInterfaceUnitTests,
    OpenTrustRegionSystemTests,
    CIntegrationTests,
    PyIntegrationTests,
)

# try to load extension module tests if available
try:
    from pyopentrustregion.extensions.common.tests import (
        CommonTests,
        CommonCInterfaceTests,
        CommonPyInterfaceTests,
    )
except AttributeError:
    pass

try:
    from pyopentrustregion.extensions.mo.tests import (
        MOTests,
        MOCInterfaceTests,
        MOPyInterfaceTests,
        MOCIntegrationTests,
        MOPyIntegrationTests,
    )
except AttributeError:
    pass

try:
    from pyopentrustregion.extensions.oao.tests import (
        OAOTests,
        OAOCInterfaceTests,
        OAOPyInterfaceTests,
        OAOCIntegrationTests,
        OAOPyIntegrationTests,
    )
except AttributeError:
    pass

try:
    from pyopentrustregion.extensions.rotation_qn.tests import (
        RotationQNTests,
        RotationQNCInterfaceTests,
        RotationQNPyInterfaceTests,
        RotationQNCIntegrationTests,
        RotationQNPyIntegrationTests,
    )
except AttributeError:
    pass

try:
    from pyopentrustregion.extensions.density_qn.tests import (
        DensityQNTests,
        DensityQNCInterfaceTests,
        DensityQNPyInterfaceTests,
        DensityQNCIntegrationTests,
        DensityQNPyIntegrationTests,
    )
except AttributeError:
    pass

try:
    from pyopentrustregion.extensions.s_gek.tests import (
        SGEKTests,
        SGEKCInterfaceTests,
        SGEKPyInterfaceTests,
        SGEKCIntegrationTests,
        SGEKPyIntegrationTests,
    )
except AttributeError:
    pass


if __name__ == "__main__":
    unittest.main(verbosity=0)
