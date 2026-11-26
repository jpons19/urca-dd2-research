"""Change only the static OPE propagator denominator, not its NN vertex.

The rho, Landau and SRC terms remain in the original Bottaro model.
This is a fixed-DD2 sensitivity test, NOT Voskresensky's collective MMU.
"""
import numpy as np


def pion_ratio(q2, m_pi, screening_ratio=1.0):
    """Return q²/(q² + m_scr²); m_scr = screening_ratio * m_pi."""
    if not np.isfinite(screening_ratio) or screening_ratio <= 0:
        raise ValueError('screening_ratio must be finite and strictly positive')
    return q2 / (q2 + (screening_ratio * m_pi)**2)


def couplings(q2, g, screening_ratio, constants):
    """Bottaro g'_q,h'_q with only the pion denominator varied.

    constants is the existing compare_compose_four module. The factor
    (f_piNN/m_pi,vac)^2 is kept fixed; so are rho and Landau terms.
    """
    rho = constants.CRHO * q2 / (q2 + constants.MRHO**2)
    pion = pion_ratio(q2, constants.MPI, screening_ratio)
    scale = (constants.FPI / constants.MPI)**2
    return g - scale*rho, -scale*(pion-rho)
