# automatic load of things related to all projects go here
# for now, using snowcrab_startup.jl to be clear

current_directory =  @__DIR__() 
print( "Current directory is: ", current_directory )

@eval using Revise
@eval using OhMyREPL
colorscheme!("Monokai24bit")

showall(x) = show(stdout, "text/plain", x)




#=

  project_directory = @__DIR__() #  same folder as the file
  push!(LOAD_PATH, project_directory)  # add the directory to the load path, so it can be found

  import Pkg  # or using Pkg

  Pkg.activate(project_directory)  # so now you activate the package
  # Pkg.activate(@__DIR__()) #  same folder as the file itself.

  Base.active_project()  # to make sure it's the package you meant to activate, print the path to console so you get a visual confirmation it's the package you meant to use

  cd( project_directory )

  print( project_directory )

=#
