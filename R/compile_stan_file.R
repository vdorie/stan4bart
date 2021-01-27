compile_stan_file <- function(file) {
  file <- sub("\\.cc$", ".stan", file)
  cppcode <- rstan::stanc(file, allow_undefined = TRUE,
                          obfuscate_model_name = FALSE)$cppcode
  cppcode <- sub("(class[[:space:]][A-Za-z_][A-Za-z0-9_]*[[:space:]])",
                 paste("#include <meta_header.hpp>\n", "\\1"), cppcode)
  
  cat("#ifndef MODELS_HPP", "#define MODELS_HPP", "#define STAN__SERVICES__COMMAND_HPP",
      "#include <rstan/rstaninc.hpp>",
      cppcode, "#endif", file = sub("\\.stan$", ".hpp", file),
      sep = "\n", append = FALSE)
  
  
  return(invisible(NULL))
}

