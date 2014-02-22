{ haskellPackages ? (import <nixpkgs> {}).haskellPackages , pkgs ? (import <nixpkgs> {})
}:

haskellPackages.cabal.mkDerivation (self: with haskellPackages; {

  pname = "MultiTimeSeries";
  version = "0.1.0.0";
  src = ./.;
  buildTools = [ cabalInstall_1_18_0_2 ];
  buildDepends = [ foldl hmatrix statistics vector ];
  doCheck = false;
  strictConfigurePhase = false;
  noHaddock = true;
  meta = {
    description = "A library for multidimensional time series analysis";
    license = self.stdenv.lib.licenses.bsd3;
    platforms = self.ghc.meta.platforms;
  };
})
