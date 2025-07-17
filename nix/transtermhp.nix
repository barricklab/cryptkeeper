{ stdenv, fetchzip  }:

stdenv.mkDerivation {
  pname = "transtermhp";
  version = "2.09";

  src = fetchzip {
    url = "https://transterm.cbcb.umd.edu/transterm_hp_v2.09.zip";
    sha256 = "sha256-Zb11M15xKjNorN0A8/wG9Yg8iWKSRgSVlCqIk1XsZts=";
  };
  patchPhase = ''
    substituteInPlace seq.cc \
      --replace "return false;" "return nullptr;"
  '';

  installPhase = ''
    mkdir -p $out/bin
    cp transterm $out/bin/
  '';
}
