using Documenter, JumpProcesses

include("pages.jl")

mathengine = MathJax3(Dict(:loader => Dict("load" => ["[tex]/require", "[tex]/mathtools"]),
                           :tex => Dict("inlineMath" => [["\$", "\$"], ["\\(", "\\)"]],
                                        "packages" => [
                                            "base",
                                            "ams",
                                            "autoload",
                                            "mathtools",
                                            "require",
                                        ])))

makedocs(sitename = "JumpProcesses.jl",
         authors = "Chris Rackauckas",
         modules = [JumpProcesses],
         clean = true,
         doctest = false,
         format = Documenter.HTML(; analytics = "UA-90474609-3",
                                  assets = ["assets/favicon.ico"],
                                  canonical = "https://jump.sciml.ai/stable/",
                                  prettyurls = (get(ENV, "CI", nothing) == "true"),
                                  mathengine),
         pages = pages)

deploydocs(repo = "github.com/SciML/JumpProcesses.jl.git";
           push_preview = true)
