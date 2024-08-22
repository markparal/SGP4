import numpy as np
from dataclasses import dataclass
import math

@dataclass
class elsetrec: # TODO Comment all variables
    satnum: str
    epochyr: int
    epochtynumrev: int
    error: int
    operationmode: str
    init: str
    method: str
    isimp: int
    aycof: float
    con41: float
    cc1: float
    cc4: float
    cc5: float
    d2: float
    d3: float
    d4: float
    delmo: float
    eta: float
    argpdot: float
    omgcof: float
    sinmao: float
    t: float
    t2cof: float
    t3cof: float
    t4cof: float
    t5cof: float
    x1mth2: float
    x7thm1: float
    mdot: float
    nodedot: float
    xlcof: float
    xmcof: float
    nodecf: float
    irez: int
    d2201: float
    d2211: float
    d3210: float
    d3222: float
    d4410: float
    d4422: float
    d5220: float
    d5232: float
    d5421: float
    d5433: float
    dedt: float
    del1: float
    del2: float
    del3: float
    didt: float
    dmdt: float
    dnodt: float
    domdt: float
    e3: float
    ee2: float
    peo: float
    pgho: float
    pho: float
    pinco: float
    plo: float
    se2: float
    se3: float
    sgh2: float
    sgh3: float
    sgh4: float
    sh2: float
    sh3: float
    si2: float
    si3: float
    sl2: float
    sl3: float
    sl4: float
    gsto: float
    xfact: float
    xgh2: float
    xgh3: float
    xgh4: float
    xh2: float
    xh3: float
    xi2: float
    xi3: float
    xl2: float
    xl3: float
    xl4: float
    xlamo: float
    zmol: float
    zmos: float
    atime: float
    xli: float
    xni: float
    a: float
    altp: float
    alta: float
    epochdays: float
    jdsatepoch: float
    jdsatepochF: float
    nddot: float
    ndot: float
    bstar: float
    rcse: float
    inclo: float
    nodeo: float
    ecco: float
    argpo: float
    mo: float
    no_kozai: float
    classification: str
    intldesg: str
    ephtype: int
    elnum: int
    revnum: int
    no_unkozai: float
    am: float
    em: float
    im: float
    Om: float
    om: float
    mm: float
    nm: float
    tumin: float
    mus: float
    radiusearthkm: float
    xke: float
    j2: float
    j3: float
    j4: float
    j3oj2: float
    dia_mm: int
    period_sec: float
    active: str
    not_orbital: str
    rcs_m2: float

class SGP4:
    def __init__(self):
        """
        TODO
        """

    def dspace(self,satrec:elsetrec,tc:float,em:float,argpm:float,inclm:float,mm:float,nodem:float,dndt:float,nm:float):
        """
        TODO Comment
        """
        fasx2 = 0.13130908
        fasx4 = 2.8843198
        fasx6 = 0.37448087
        g22 = 5.7686396
        g32 = 0.95240898
        g44 = 1.8014998
        g52 = 1.0508330
        g54 = 4.4108898
        rptim = 4.37526908801129966e-3
        stepp = 720.0
        stepn = -720.0
        step2 = 259200.0

        # Calculate deep space resonance effects
        dndt = 0.0
        theta = np.fmod(satrec.gsto+tc*rptim,2*np.pi)
        em = em + satrec.dedt * satrec.t

        inclm = inclm + satrec.didt * satrec.t
        argpm = argpm + satrec.domdt * satrec.t
        nodem = nodem + satrec.dnodt * satrec.t
        mm = mm + satrec.dmdt * satrec.t

        # sgp4fix take out satrec.atime = 0.0 and fix for faster operation
        ft = 0.0
        if (satrec.irez != 0):
            # sgp4fix streamline check
            if ((satrec.atime == 0.0) or (satrec.t * satrec.atime <= 0.0) or (np.abs(satrec.t) < np.abs(satrec.atime))):
                satrec.atime = 0.0
                satrec.xni = satrec.no_unkozai
                satrec.xli = satrec.xlamo

            # sgp4fix move check outside loop
            if (satrec.t > 0.0):
                delt = stepp
            else:
                delt = stepn

            iretn = 381 # added for do loop
            iret = 0 # added for loop
            while (iretn == 381):
                # near - synchronous resonance terms
                if (satrec.irez != 2):
                    xndt = satrec.del1 * np.sin(satrec.xli - fasx2) + satrec.del2 * np.sin(2.0 * (satrec.xli - fasx4)) + satrec.del3 * np.sin(3.0 * (satrec.xli - fasx6))
                    xldot = satrec.xni + satrec.xfact
                    xnddt = satrec.del1 * np.cos(satrec.xli - fasx2) + 2.0 * satrec.del2 * np.cos(2.0 * (satrec.xli - fasx4)) + 3.0 * satrec.del3 * np.cos(3.0 * (satrec.xli - fasx6))
                    xnddt = xnddt * xldot
                else:
                    # near - half-day resonance terms
                    xomi = satrec.argpo + satrec.argpdot * satrec.atime
                    x2omi = xomi + xomi
                    x2li = satrec.xli + satrec.xli
                    xndt = satrec.d2201 * np.sin(x2omi + satrec.xli - g22) + satrec.d2211 * np.sin(satrec.xli - g22) + \
                        satrec.d3210 * np.sin(xomi + satrec.xli - g32) + satrec.d3222 * np.sin(-xomi + satrec.xli - g32) + \
                        satrec.d4410 * np.sin(x2omi + x2li - g44) + satrec.d4422 * np.sin(x2li - g44) + \
                        satrec.d5220 * np.sin(xomi + satrec.xli - g52) + satrec.d5232 * np.sin(-xomi + satrec.xli - g52) + \
                        satrec.d5421 * np.sin(xomi + x2li - g54) + satrec.d5433 * np.sin(-xomi + x2li - g54)
                    xldot = satrec.xni + satrec.xfact
                    xnddt = satrec.d2201 * np.cos(x2omi + satrec.xli - g22) + satrec.d2211 * np.cos(satrec.xli - g22) + \
                        satrec.d3210 * np.cos(xomi + satrec.xli - g32) + satrec.d3222 * np.cos(-xomi + satrec.xli - g32) + \
                        satrec.d5220 * np.cos(xomi + satrec.xli - g52) + satrec.d5232 * np.cos(-xomi + satrec.xli - g52) + \
                        2.0 * (satrec.d4410 * np.cos(x2omi + x2li - g44) + \
                        satrec.d4422 * np.cos(x2li - g44) + satrec.d5421 * np.cos(xomi + x2li - g54) + \
                        satrec.d5433 * np.cos(-xomi + x2li - g54))
                    xnddt = xnddt * xldot

                # integrator
                # sgp4fix move end checks to end of routine
                if (np.abs(satrec.t - satrec.atime) >= stepp):
                    iret = 0
                    iretn = 381
                else:
                    ft = satrec.t - satrec.atime
                    iretn = 0

                if (iretn == 381):
                    xli = xli + xldot * delt + xndt * step2
                    satrec.xni = satrec.xni + xndt * delt + xnddt * step2
                    satrec.atime = satrec.atime + delt

            nm = satrec.xni + xndt * ft + xnddt * ft * ft * 0.5
            xl = xli + xldot * ft + xndt * ft * ft * 0.5
            if (satrec.irez != 1):
                mm = xl - 2.0 * nodem + 2.0 * theta
                dndt = nm - satrec.no_unkozai
            else:
                mm = xl - nodem - argpm + theta
                dndt = nm - satrec.no_unkozai
            nm = satrec.no_unkozai + dndt
        
        return satrec,em,argpm,inclm,mm,nodem,dndt,nm
    
    def dpper(satrec:elsetrec,init:str,ep:float,inclp:float,nodep:float,argpp:float,mp:float):
        """
        """

        # constants
        zns = 1.19459e-5
        zes = 0.01675
        znl = 1.5835218e-4
        zel = 0.05490

		# calculate time varying periodics
        zm = satrec.zmos + zns * satrec.t
        # be sure that the initial call has time set to zero
        if (init == 'y'):
            zm = satrec.zmos
        zf = zm + 2.0 * zes * np.sin(zm)
        sinzf = np.sin(zf)
        f2 = 0.5 * sinzf * sinzf - 0.25
        f3 = -0.5 * sinzf * np.cos(zf)
        ses = satrec.se2* f2 + satrec.se3 * f3
        sis = satrec.si2 * f2 + satrec.si3 * f3
        sls = satrec.sl2 * f2 + satrec.sl3 * f3 + satrec.sl4 * sinzf
        sghs = satrec.sgh2 * f2 + satrec.sgh3 * f3 + satrec.sgh4 * sinzf
        shs = satrec.sh2 * f2 + satrec.sh3 * f3
        zm = satrec.zmol + znl * satrec.t
        if (init == 'y'):
            zm = satrec.zmol
        zf = zm + 2.0 * zel * np.sin(zm)
        sinzf = np.sin(zf)
        f2 = 0.5 * sinzf * sinzf - 0.25
        f3 = -0.5 * sinzf * np.cos(zf)
        sel = satrec.ee2 * f2 + satrec.e3 * f3
        sil = satrec.xi2 * f2 + satrec.xi3 * f3
        sll = satrec.xl2 * f2 + satrec.xl3 * f3 +satrec.xl4 * sinzf
        sghl = satrec.xgh2 * f2 + satrec.xgh3 * f3 + satrec.xgh4 * sinzf
        shll = satrec.xh2 * f2 + satrec.xh3 * f3
        pe = ses + sel
        pinc = sis + sil
        pl = sls + sll
        pgh = sghs + sghl
        ph = shs + shll

        if (init == 'n'):
            pe = pe - satrec.peo
            pinc = pinc - satrec.pinco
            pl = pl - satrec.plo
            pgh = pgh - satrec.pgho
            ph = ph - satrec.pho
            inclp = inclp + pinc
            ep = ep + pe
            sinip = np.sin(inclp)
            cosip = np.cos(inclp)

            # apply periodics directly
			#  sgp4fix for lyddane choice
			#  strn3 used original inclination - this is technically feasible
			#  gsfc used perturbed inclination - also technically feasible
			#  probably best to readjust the 0.2 limit value and limit discontinuity
			#  0.2 rad = 11.45916 deg
			#  use next line for original strn3 approach and original inclination
			#  if (inclo >= 0.2)
			#  use next line for gsfc version and perturbed inclination
            if (inclp >= 0.2):
                ph = ph / sinip
                pgh = pgh - cosip * ph
                argpp = argpp + pgh
                nodep = nodep + ph
                mp = mp + pl
            else:
				# apply periodics with lyddane modification
                sinop = np.sin(nodep)
                cosop = np.cos(nodep)
                alfdp = sinip * sinop
                betdp = sinip * cosop
                dalf = ph * cosop + pinc * cosip * sinop
                dbet = -ph * sinop + pinc * cosip * cosop
                alfdp = alfdp + dalf
                betdp = betdp + dbet
                nodep = np.fmod(nodep, 2*np.pi)
                #  sgp4fix for afspc written intrinsic functions
                # nodep used without a trigonometric function ahead
                if ((nodep < 0.0) and (satrec.operationmode == 'a')):
                    nodep = nodep + 2*np.pi
                xls = mp + argpp + cosip * nodep
                dls = pl + pgh - pinc * nodep * sinip
                xls = xls + dls
                xnoh = nodep
                nodep = math.atan2(alfdp, betdp)
                #  sgp4fix for afspc written intrinsic functions
                # nodep used without a trigonometric function ahead
                if ((nodep < 0.0) and (satrec.operationmode == 'a')):
                    nodep = nodep + 2*np.pi
                if (np.fabs(xnoh - nodep) > np.pi):
                    if (nodep < xnoh):
                        nodep = nodep + 2*np.pi
                    else:
                        nodep = nodep - 2*np.pi
                mp = mp + pl
                argpp = xls - mp - cosip * nodep

        return satrec,ep,inclp,nodep,argpp,mp
    
    def sgp4_prop(self,satrec:elsetrec,tsince:float,r:np.ndarray,v:np.ndarray):
        """
        This procedure is the sgp4 prediction model from space command. this is an
        updated and combined version of sgp4 and sdp4, which were originally published 
        separately in spacetrack report #3. this version follows the methodology from 
        the aiaa paper (2006) describing the history and development of the code.

        Original Author: David Vallado
        Python Port Author: Mark Paral

        Parameters
        ----------
        satrec: elsetrec
            TODO
        tsince: float
            TODO
        r: ndarray
            TODO
        v: ndarray
            TODO

        Returns
        -------
        TODO
        """

        # Set mathematical constants TODO description
        temp4 = 1.5e-12
        x2o3 = 2.0 / 3.0

        vkmpersec = satrec.radiusearthkm * satrec.xke / 60.0

        satrec.t = tsince
        satrec.error = 0

        xmdf = satrec.mo + satrec.mdot * satrec.t
        argpdf = satrec.argpo + satrec.argpdot * satrec.t
        nodedf = satrec.nodeo + satrec.nodedot * satrec.t
        argpm = argpdf
        mm = xmdf
        t2 = satrec.t * satrec.t
        nodem = nodedf + satrec.nodecf * t2
        tempa = 1.0 - satrec.cc1 * satrec.t
        tempe = satrec.bstar * satrec.cc4 * satrec.t
        templ = satrec.t2cof * t2

        if (satrec.isimp != 1):
            delomg = satrec.omgcof * satrec.t
            delmtemp = 1.0 + satrec.eta * np.cos(xmdf)
            delm = satrec.xmcof * (delmtemp * delmtemp * delmtemp - satrec.delmo)
            temp = delomg + delm
            mm = xmdf + temp
            argpm = argpdf - temp
            t3 = t2 * satrec.t
            t4 = t3 * satrec.t
            tempa = tempa - satrec.d2 * t2 - satrec.d3 * t3 - satrec.d4 * t4
            tempe = tempe + satrec.bstar * satrec.cc5 * (np.sin(mm) - satrec.sinmao)
            templ = templ + satrec.t3cof * t3 + t4 * (satrec.t4cof + satrec.t * satrec.t5cof)

        nm = satrec.no_unkozai
        em = satrec.ecco
        inclm = satrec.inclo

        if (satrec.method == 'd'):
            tc = satrec.t
            dndt = 0.
            satrec,em,argpm,inclm,mm,nodem,dndt,nm = self.dspace(self,satrec,tc,em,argpm,inclm,mm,nodem,dndt,nm)

        if (nm <= 0.0):
            satrec.error = 2
            # sgp4fix add return
            return False
        
        am = pow((satrec.xke / nm), x2o3) * tempa * tempa
        nm = satrec.xke / pow(am, 1.5)
        em = em - tempe

        if ((em >= 1.0) or (em < -0.001)):
            satrec.error = 1
            # sgp4fix to return if there is an error in eccentricity
            return False
        
        # sgp4fix fix tolerance to avoid a divide by zero
        if (em < 1.0e-6):
            em = 1.0e-6
        mm = mm + satrec.no_unkozai * templ
        xlm = mm + argpm + nodem
        emsq = em * em
        temp = 1.0 - emsq

        nodem = np.fmod(nodem, 2*np.pi)
        argpm = np.fmod(argpm, 2*np.pi)
        xlm = np.fmod(xlm, 2*np.pi)
        mm = np.fmod(xlm - argpm - nodem, 2*np.pi)

        # sgp4fix recover singly averaged mean elements
        satrec.am = am
        satrec.em = em
        satrec.im = inclm
        satrec.Om = nodem
        satrec.om = argpm
        satrec.mm = mm
        satrec.nm = nm

        # compute extra mean quantities
        sinim = np.sin(inclm)
        cosim = np.cos(inclm)

		# add lunar-solar periodics
        ep = em
        xincp = inclm
        argpp = argpm
        nodep = nodem
        mp = mm
        sinip = sinim
        cosip = cosim
        if (satrec.method == 'd'):
            init = 'n'
            satrec,ep,inclp,nodep,argpp,mp = self.dpper(satrec,init,ep,inclp,nodep,argpp,mp)
            if (xincp < 0.0):
                xincp = -xincp
                nodep = nodep + np.pi
                argpp = argpp - np.pi
            
            if ((ep < 0.0) or (ep > 1.0)):
                satrec.error = 3
				# sgp4fix add return
                return False

		# long period periodics 
        if (satrec.method == 'd'):
            sinip = np.sin(xincp)
            cosip = np.cos(xincp)
            satrec.aycof = -0.5*satrec.j3oj2*sinip
            # sgp4fix for divide by zero for xincp = 180 deg
            if (np.fabs(cosip + 1.0) > 1.5e-12):
                satrec.xlcof = -0.25 * satrec.j3oj2 * sinip * (3.0 + 5.0 * cosip) / (1.0 + cosip)
            else:
                satrec.xlcof = -0.25 * satrec.j3oj2 * sinip * (3.0 + 5.0 * cosip) / temp4

        axnl = ep * np.cos(argpp)
        temp = 1.0 / (am * (1.0 - ep * ep))
        aynl = ep* np.sin(argpp) + temp * satrec.aycof
        xl = mp + argpp + nodep + temp * satrec.xlcof * axnl

        # solve kepler's equation
        u = np.fmod(xl - nodep, 2*np.pi)
        eo1 = u
        tem5 = 9999.9
        ktr = 1
        #   sgp4fix for kepler iteration
        #   the following iteration needs better limits on corrections
        while ((np.fabs(tem5) >= 1.0e-12) and (ktr <= 10)):
            sineo1 = np.sin(eo1)
            coseo1 = np.cos(eo1)
            tem5 = 1.0 - coseo1 * axnl - sineo1 * aynl
            tem5 = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5
            if (np.fabs(tem5) >= 0.95):
                tem5 = 0.95 if tem5 > 0.0 else -0.95
            eo1 = eo1 + tem5
            ktr = ktr + 1

        # short period preliminary quantities
        ecose = axnl*coseo1 + aynl*sineo1
        esine = axnl*sineo1 - aynl*coseo1
        el2 = axnl*axnl + aynl*aynl
        pl = am*(1.0 - el2)
        if (pl < 0.0):
            satrec.error = 4
            # sgp4fix add return
            return False
        else:
            rl = am * (1.0 - ecose)
            rdotl = np.sqrt(am) * esine / rl
            rvdotl = np.sqrt(pl) / rl
            betal = np.sqrt(1.0 - el2)
            temp = esine / (1.0 + betal)
            sinu = am / rl * (sineo1 - aynl - axnl * temp)
            cosu = am / rl * (coseo1 - axnl + aynl * temp)
            su = math.atan2(sinu, cosu)
            sin2u = (cosu + cosu) * sinu
            cos2u = 1.0 - 2.0 * sinu * sinu
            temp = 1.0 / pl
            temp1 = 0.5 * satrec.j2 * temp
            temp2 = temp1 * temp

            # update for short period periodics
            if (satrec.method == 'd'):
                cosisq = cosip * cosip
                satrec.con41 = 3.0*cosisq - 1.0
                satrec.x1mth2 = 1.0 - cosisq
                satrec.x7thm1 = 7.0*cosisq - 1.0
            mrt = rl * (1.0 - 1.5 * temp2 * betal * satrec.con41) + 0.5 * temp1 * satrec.x1mth2 * cos2u
            su = su - 0.25 * temp2 * satrec.x7thm1 * sin2u
            xnode = nodep + 1.5 * temp2 * cosip * sin2u
            xinc = xincp + 1.5 * temp2 * cosip * sinip * cos2u
            mvt = rdotl - nm * temp1 * satrec.x1mth2 * sin2u / satrec.xke
            rvdot = rvdotl + nm * temp1 * (satrec.x1mth2 * cos2u + 1.5 * satrec.con41) / satrec.xke

            # orientation vectors
            sinsu = np.sin(su)
            cossu = np.cos(su)
            snod = np.sin(xnode)
            cnod = np.cos(xnode)
            sini = np.sin(xinc)
            cosi = np.cos(xinc)
            xmx = -snod * cosi
            xmy = cnod * cosi
            ux = xmx * sinsu + cnod * cossu
            uy = xmy * sinsu + snod * cossu
            uz = sini * sinsu
            vx = xmx * cossu - cnod * sinsu
            vy = xmy * cossu - snod * sinsu
            vz = sini * cossu

            # position and velocity (in km and km/sec)
            r[0] = (mrt * ux)* satrec.radiusearthkm
            r[1] = (mrt * uy)* satrec.radiusearthkm
            r[2] = (mrt * uz)* satrec.radiusearthkm
            v[0] = (mvt * ux + rvdot * vx) * vkmpersec
            v[1] = (mvt * uy + rvdot * vy) * vkmpersec
            v[2] = (mvt * uz + rvdot * vz) * vkmpersec

		# sgp4fix for decaying satellites
        if (mrt < 1.0):
            satrec.error = 6
            return False
        
        return True