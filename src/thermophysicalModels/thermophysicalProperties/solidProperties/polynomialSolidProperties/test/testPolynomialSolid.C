/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    testPolynomialSolid

Description
    Unit test for polynomialSolid thermophysical properties library.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polynomialSolid.H"
#include "IStringStream.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList args(argc, argv);

    Info<< "Testing polynomialSolid implementation..." << nl << endl;

    IStringStream is
    (
        "polynomialSolid\n"
        "{\n"
        "    Hf          0.0;\n"
        "    emissivity  0.9;\n"
        "    rhoCoeffs    (2500 -0.1 0 0 0 0 0 0);\n"
        "    CpCoeffs     (800 0.5 0 0 0 0 0 0);\n"
        "    kappaCoeffs  (1.5 0 0 0 0 0 0 0);\n"
        "}\n"
    );

    dictionary dict(is);
    const dictionary& solidDict = dict.subDict("polynomialSolid");

    // 1. Runtime selection test
    autoPtr<solidProperties> solid = solidProperties::New(solidDict);

    if (!solid.valid())
    {
        FatalErrorInFunction << "Failed to construct solidProperties" << exit(FatalError);
    }
    Info<< "1. Runtime Selection:" << nl
        << "   Type: " << solid->type() << " (OK)" << nl << endl;

    const polynomialSolid& ps = refCast<const polynomialSolid>(solid());

    // 2. Direct property tests
    Info<< "2. Evaluating Polynomial Values:" << nl;

    // rho(T) = 2500 - 0.1*T
    // At T=300: 2500 - 30 = 2470
    // At T=500: 2500 - 50 = 2450
    const scalar rho300 = ps.rho(300);
    const scalar rho500 = ps.rho(500);
    Info<< "   rho(300) = " << rho300 << " (expected 2470)" << nl
        << "   rho(500) = " << rho500 << " (expected 2450)" << nl;

    if (mag(rho300 - 2470) > 1e-6 || mag(rho500 - 2450) > 1e-6)
    {
        FatalErrorInFunction << "rho evaluation mismatch!" << exit(FatalError);
    }

    // Cp(T) = 800 + 0.5*T
    // At T=300: 800 + 150 = 950
    // At T=500: 800 + 250 = 1050
    const scalar Cp300 = ps.Cp(300);
    const scalar Cp500 = ps.Cp(500);
    Info<< "   Cp(300)  = " << Cp300 << " (expected 950)" << nl
        << "   Cp(500)  = " << Cp500 << " (expected 1050)" << nl;

    if (mag(Cp300 - 950) > 1e-6 || mag(Cp500 - 1050) > 1e-6)
    {
        FatalErrorInFunction << "Cp evaluation mismatch!" << exit(FatalError);
    }

    const scalar kappa300 = ps.kappa(300);
    Info<< "   kappa(300) = " << kappa300 << " (expected 1.5)" << nl;
    if (mag(kappa300 - 1.5) > 1e-6)
    {
        FatalErrorInFunction << "kappa evaluation mismatch!" << exit(FatalError);
    }

    // 3. Enthalpy integration test: hs(T) = int_{Tstd}^T (800 + 0.5*t) dt
    const scalar hs500 = ps.hs(500);
    const scalar expectedHs500 = 800.0*(500.0 - Tstd) + 0.25*(sqr(500.0) - sqr(Tstd));
    Info<< "   hs(500)  = " << hs500 << " (expected " << expectedHs500 << ") J/kg" << nl << endl;
    if (mag(hs500 - expectedHs500) > 1e-4)
    {
        FatalErrorInFunction << "hs evaluation mismatch!" << exit(FatalError);
    }

    // 4. Default 0-argument overloads (evaluating at Tstd = 298.15 K)
    Info<< "3. Base Class 0-Argument Overloads (at Tstd = " << Tstd << " K):" << nl;
    Info<< "   rho()    = " << ps.rho() << " (rho(Tstd) = " << ps.rho(Tstd) << ")" << nl
        << "   Cp()     = " << ps.Cp() << " (Cp(Tstd) = " << ps.Cp(Tstd) << ")" << nl
        << "   kappa()  = " << ps.kappa() << " (kappa(Tstd) = " << ps.kappa(Tstd) << ")" << nl;

    if (mag(ps.rho() - ps.rho(Tstd)) > 1e-6 || mag(ps.Cp() - ps.Cp(Tstd)) > 1e-6)
    {
        FatalErrorInFunction << "0-argument overload mismatch!" << exit(FatalError);
    }
    Info<< "   Overloads match Tstd evaluation (OK)" << nl << endl;

    // 5. Clone test (validates deep copy)
    Info<< "4. Clone and Copy Integrity:" << nl;
    autoPtr<solidProperties> cloned = solid->clone();
    const polynomialSolid& psCloned = refCast<const polynomialSolid>(cloned());

    Info<< "   Cloned solid rho(500): " << psCloned.rho(500) << nl
        << "   Original solid rho(500) after clone: " << ps.rho(500) << nl;

    if (mag(psCloned.rho(500) - 2450) > 1e-6 || mag(ps.rho(500) - 2450) > 1e-6)
    {
        FatalErrorInFunction << "Clone check failed!" << exit(FatalError);
    }
    Info<< "   Copy verified (OK)" << nl << endl;

    // 6. Serialization test
    Info<< "5. Dictionary Serialization (write):" << nl;
    ps.write(Info);

    Info<< nl << "All tests completed successfully!" << endl;
    return 0;
}
