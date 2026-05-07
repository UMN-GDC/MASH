require("dap-python").setup("/usr/bin/python")

table.insert(require("dap").configurations.python, {
    type = "python",
    request = "launch",
    name = "Module",
    console = "integratedTerminal",
    module = "src", -- edit this to be your app's main module
    cwd = "${workspaceFolder}",
})
