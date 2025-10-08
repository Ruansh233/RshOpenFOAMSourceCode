    Info<< nl << "Create Times for case: " << masterCase.name() << endl;

    const fileName masterCasePath = masterCase.path();
    const fileName masterCaseName = masterCase.name();

    Time runTimeMaster
    (
        Time::controlDictName,
        masterCasePath,
        masterCaseName
    );
    runTimeMaster.functionObjects().off();
