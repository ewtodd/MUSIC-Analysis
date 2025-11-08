{
  description = "ROOT Waveform Analysis Framework";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    toolkit.url = "github:ewtodd/Nuclear-Measurement-Toolkit";
  };

  outputs = { self, nixpkgs, flake-utils, toolkit }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
        nm-toolkit = toolkit.packages.${system}.default;
      in {
        devShells.default = pkgs.mkShell {
          buildInputs = with pkgs; [ root gnumake pkg-config clang-tools ];

          shellHook = "";
        };
      });
}
