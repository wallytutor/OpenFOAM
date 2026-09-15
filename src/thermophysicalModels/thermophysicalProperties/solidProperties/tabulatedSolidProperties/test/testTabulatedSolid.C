/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    testTabulatedSolid

Description
    Unit test for tabulatedSolid thermophysical properties library.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "tabulatedSolid.H"
#include "IStringStream.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList args(argc, argv);

    Info<< "Testing tabulatedSolid implementation..." << nl << endl;

    IStringStream is
    (
        "tabulatedSolid\n"
        "{\n"
        "    Hf          0.0;\n"
        "    emissivity  0.9;\n"
        "    rho table\n"
        "    (\n"
        "        (300  2500)\n"
        "        (600  2450)\n"
        "        (1200 2380)\n"
        "    );\n"
        "    Cp table\n"
        "    (\n"
        "        (300  800)\n"
        "        (500  950)\n"
        "        (800  1100)\n"
        "        (1200 1250)\n"
        "    );\n"
        "    kappa table\n"
        "    (\n"
        "        (300  1.5)\n"
        "        (500  1.5)\n"
        "        (800  1.5)\n"
        "        (1200 1.5)\n"
        "    );\n"
        "}\n"
    );

    dictionary dict(is);
    const dictionary& solidDict = dict.subDict("tabulatedSolid");

    // 1. Runtime selection test
    autoPtr<solidProperties> solid = solidProperties::New(solidDict);

    if (!solid.valid())
    {
        FatalErrorInFunction << "Failed to construct solidProperties" << exit(FatalError);
    }
    Info<< "1. Runtime Selection:" << nl
        << "   Type: " << solid->type() << " (OK)" << nl << endl;

    const tabulatedSolid& ts = refCast<const tabulatedSolid>(solid());

    // 2. Direct property tests
    Info<< "2. Evaluating Tabulated & Interpolated Values:" << nl;

    const scalar rho300 = ts.rho(300);
    const scalar rho450 = ts.rho(450); // Linear interpolation: (2500 + 2450)/2 = 2475
    const scalar rho600 = ts.rho(600);
    Info<< "   rho(300) = " << rho300 << " (expected 2500)" << nl
        << "   rho(450) = " << rho450 << " (expected 2475)" << nl
        << "   rho(600) = " << rho600 << " (expected 2450)" << nl;

    if (mag(rho300 - 2500) > 1e-6 || mag(rho450 - 2475) > 1e-6 || mag(rho600 - 2450) > 1e-6)
    {
        FatalErrorInFunction << "rho evaluation mismatch!" << exit(FatalError);
    }

    const scalar Cp300 = ts.Cp(300);
    const scalar Cp500 = ts.Cp(500);
    Info<< "   Cp(300)  = " << Cp300 << " (expected 800)" << nl
        << "   Cp(500)  = " << Cp500 << " (expected 950)" << nl;

    if (mag(Cp300 - 800) > 1e-6 || mag(Cp500 - 950) > 1e-6)
    {
        FatalErrorInFunction << "Cp evaluation mismatch!" << exit(FatalError);
    }

    const scalar kappa300 = ts.kappa(300);
    Info<< "   kappa(300) = " << kappa300 << " (expected 1.5)" << nl;
    if (mag(kappa300 - 1.5) > 1e-6)
    {
        FatalErrorInFunction << "kappa evaluation mismatch!" << exit(FatalError);
    }

    // 3. Enthalpy integration test
    const scalar hs500 = ts.hs(500);
    Info<< "   hs(500)  = " << hs500 << " J/kg" << nl << endl;

    // 4. Default 0-argument overloads (evaluating at Tstd = 298.15 K)
    Info<< "3. Base Class 0-Argument Overloads (at Tstd = " << Tstd << " K):" << nl;
    Info<< "   rho()    = " << ts.rho() << " (rho(Tstd) = " << ts.rho(Tstd) << ")" << nl
        << "   Cp()     = " << ts.Cp() << " (Cp(Tstd) = " << ts.Cp(Tstd) << ")" << nl
        << "   kappa()  = " << ts.kappa() << " (kappa(Tstd) = " << ts.kappa(Tstd) << ")" << nl;

    if (mag(ts.rho() - ts.rho(Tstd)) > 1e-6 || mag(ts.Cp() - ts.Cp(Tstd)) > 1e-6)
    {
        FatalErrorInFunction << "0-argument overload mismatch!" << exit(FatalError);
    }
    Info<< "   Overloads match Tstd evaluation (OK)" << nl << endl;

    // 5. Clone test (validates autoPtr deep-copy)
    Info<< "4. Clone and Deep Copy Integrity:" << nl;
    autoPtr<solidProperties> cloned = solid->clone();
    const tabulatedSolid& tsCloned = refCast<const tabulatedSolid>(cloned());

    Info<< "   Cloned solid rho(600): " << tsCloned.rho(600) << nl
        << "   Original solid rho(600) after clone: " << ts.rho(600) << nl;

    if (mag(tsCloned.rho(600) - 2450) > 1e-6 || mag(ts.rho(600) - 2450) > 1e-6)
    {
        FatalErrorInFunction << "Clone check failed: autoPtr was moved instead of deep copied!"
            << exit(FatalError);
    }
    Info<< "   Deep copy verified (OK)" << nl << endl;

    // 6. Serialization test
    Info<< "5. Dictionary Serialization (write):" << nl;
    ts.write(Info);

    Info<< nl << "All tests completed successfully!" << endl;
    return 0;
}
