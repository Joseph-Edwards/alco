LoadPackage("GapDoc");

pkg := "ALCO";

dir := DirectoriesPackageLibrary(pkg, "doc")[1];

MakeGAPDocDoc(dir, "alco.xml", [], pkg);
CopyHTMLStyleFiles( "doc" );

QUIT;
