{ stdenv, fetchzip  }:

stdenv.mkDerivation {
  pname = "transtermhp";
  version = "2.09";

  src = fetchzip {
    url = "https://depot.galaxyproject.org/software/transtermhp/transtermhp_2.09_src_all.zip";
    hash = "sha256-Zb11M15xKjNorN0A8/wG9Yg8iWKSRgSVlCqIk1XsZts=";
  };

  patchPhase = ''
    substituteInPlace seq.cc \
      --replace "return false;" "return nullptr;"
  '';

  installPhase = ''
    mkdir -p $out/bin $out/data
    cp transterm $out/bin/
    cp expterm.dat $out/data/
  '';
}
