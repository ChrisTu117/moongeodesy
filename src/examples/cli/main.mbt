///|
/// Command-line demo. Run with `moon run src/examples/cli`.
fn main {
  match execute(@env.args()) {
    Ok(text) => println(text)
    Err(text) => {
      eprintln(text)
      exit_with(1)
    }
  }
}

///|
#cfg(target="wasm")
fn exit_with(code : Int32) -> Unit = "wasi_snapshot_preview1" "proc_exit"

///|
#cfg(target="wasm-gc")
fn exit_with(code : Int32) -> Unit = "__moonbit_sys_unstable" "exit"

///|
#cfg(target="js")
extern "js" fn exit_with(code : Int32) -> Unit =
  #| function(code) {
  #|   if (typeof process !== "undefined" && typeof process.exit === "function") {
  #|     process.exit(code);
  #|   }
  #| }

///|
#cfg(any(target="native", target="llvm"))
extern "c" fn exit_with(code : Int32) = "exit"
