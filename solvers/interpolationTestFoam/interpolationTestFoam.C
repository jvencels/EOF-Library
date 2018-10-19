/*---------------------------------------------------------------------------*\

Author
    Juris Vencels

Description
    Solver for testing EOF-Library interpolation and communication.

\*---------------------------------------------------------------------------*/

#include "fvOptions.H"
#include "Elmer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFvOptions.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "Reading fields\n" << endl;

    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Receive fields from Elmer
    Elmer<fvMesh> receiving(mesh,-1); // 1=send, -1=receive
    receiving.sendStatus(0); // 1=ok, 0=lastIter, -1=error
    receiving.recvScalar(T);

/*---------------------------------------------------------------------------*\
    volScalarField realT
    (
        IOobject
        (
            "realT",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    forAll(realT.mesh().cells(),cellI) 
    {
        scalar x = realT.mesh().C()[cellI].x();
        scalar y = realT.mesh().C()[cellI].y();
        scalar z = realT.mesh().C()[cellI].z();

        realT[cellI] = 2*Foam::sin(2*Foam::constant::mathematical::pi*x)*Foam::cos(2*Foam::constant::mathematical::pi*y);
    }

    fvOptions.correct(realT);
    realT.write();

    Info<< "Elmer2OpenFOAM error: " << Foam::sqrt(gSum((sqr(realT-T))()))/Foam::sqrt(gSum((sqr(realT))())) << "\n" << endl;
\*---------------------------------------------------------------------------*/

    fvOptions.correct(T);
    T.write();

    // Send fields to Elmer
    Elmer<fvMesh> sending(mesh,1); // 1=send, -1=receive
    sending.sendStatus(0); // 1=ok, 0=lastIter, -1=error
    sending.sendScalar(T);
    //sending.sendScalar(realT);

    return 0;
}


// ************************************************************************* //
