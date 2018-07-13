function out = runTests( )
    out = matlab.unittest.TestSuite.fromPackage('nirs.testing.unittests').run;
end

